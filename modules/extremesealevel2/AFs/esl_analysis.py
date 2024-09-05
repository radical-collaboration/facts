'''
@author: Tim Hermans
t(dot)h(dot)j(dot)hermans@uu(dot)nl
'''
import numpy as np
import numpy.matlib
import pandas as pd
import xarray as xr
from scipy.stats import genpareto
from utils import mindist
from tqdm import tqdm
from preprocessing import extract_GESLA2_locations, extract_GESLA3_locations, ingest_GESLA2_files, ingest_GESLA3_files 
from preprocessing import detrend_gesla_dfs, deseasonalize_gesla_dfs, subtract_amean_from_gesla_dfs
from preprocessing import detrend_ds, deseasonalize_ds, subtract_amean_from_ds
from I_O import open_gtsm_waterlevels
import os

def ESL_stats_from_gtsm_dmax(gtsm_path,input_locations,preproc_settings,match_dist_limit,output_dir=None):
    
    print('Warning: Several preprocessing options applicable to GESLA data are bypassed for analyzing GTSM data, such as minimum length, MSL reference, and modeling events below the GPD threshold.')
    
    assert ((preproc_settings['extremes_threshold']>0) & (preproc_settings['extremes_threshold']<100))
    
    if preproc_settings['resample_freq'] != 'D_max':
        raise Exception('Configured resample frequency not yet implemented. Assuming daily maxima Gtsm data.')
    if preproc_settings['declus_method'] != 'rolling_max':
        raise Exception('Configured declustering method not yet implemented.')      
    
    print('Opening GTSM timeseries near input_locations sites.')
    gtsm = open_gtsm_waterlevels(gtsm_path,input_locations,match_dist_limit)
    
    print('Preprocessing GTSM timeseries & finding extremes.')
    #preproc options:
    if preproc_settings['detrend']:
        gtsm = detrend_ds(gtsm,'waterlevel')
    if preproc_settings['deseasonalize']:
        gtsm = deseasonalize_ds(gtsm,'waterlevel')
    if preproc_settings['subtract_amean']:
        gtsm = subtract_amean_from_ds(gtsm,'waterlevel')
        
    threshold= gtsm.waterlevel.quantile(preproc_settings['extremes_threshold']/100,dim='time')
    extremes = gtsm.waterlevel.where(gtsm.waterlevel>=threshold)
    extremes_declustered = extremes.where((extremes==extremes.rolling(time=preproc_settings['declus_window'],center=True,min_periods=1).max()) & (np.isfinite(extremes))) #note this only works if resolution is daily!!
    
    if preproc_settings['store_esls']:
        extremes_declustered.to_netcdf(os.path.join(output_dir,'GTSM_dmax_declustered_extremes.nc'),mode='w')
            
    avg_extr_pyear = 365.25 * np.isfinite(extremes_declustered).sum(dim='time')/len(gtsm.time)
    
    print('Fitting GPD to extremes.')
    gpd_params,nlogl,cov_mats = xr.apply_ufunc(fit_gpd_to_extremes_ds,
                   (extremes_declustered-threshold),
                   input_core_dims=[["time"]],  # list with one entry per arg
                   output_core_dims=[['params'],[],['i','j']],
                   vectorize=True,
                   dask="parallelized")

    esl_statistics = gpd_params.to_dataset(name='gpd_params')
    esl_statistics = esl_statistics.assign_coords({'params':['scale','shape']})
    esl_statistics['scale'] = esl_statistics.gpd_params.sel(params='scale')
    esl_statistics['shape'] = esl_statistics.gpd_params.sel(params='shape')
    esl_statistics['loc'] = threshold
    esl_statistics['nlogl'] = nlogl
    esl_statistics['cov'] = cov_mats
    esl_statistics['avg_extr_pyear'] = avg_extr_pyear
    esl_statistics = esl_statistics.drop(['gpd_params','params'])
    
    return esl_statistics


def ESL_stats_from_raw_GESLA(gesla_version,path_to_gesla,input_locations,preproc_settings,match_dist_limit,output_dir=None):
    ''' For input_locations sites, try to find nearest GESLA record within "maxdist" that fulfills the criteria for being included set in "cfg".'''    

    #get GESLA locations
    if gesla_version == 'gesla2':
        (station_names, station_lats, station_lons, station_filenames) = extract_GESLA2_locations(path_to_gesla,preproc_settings['min_yrs'])#extract_GESLA2_locations('/Volumes/Naamloos/PhD_Data/GESLA2/private_14032017_public_110292018/')
    elif gesla_version == 'gesla3':
        (station_names, station_lats, station_lons, station_filenames) = extract_GESLA3_locations(path_to_gesla,preproc_settings['min_yrs'])
    else:
        raise Exception('Data source not recognized.')
    
    min_idx = [mindist(x,y,station_lats,station_lons,match_dist_limit) for x,y in zip(input_locations.lat.values, input_locations.lon.values)] #sites within maximum distance from input_locations points

    matched_filenames = []
    sites_with_esl = []
    for i in np.arange(len(input_locations.locations)): #loop over input_locations sites
        if min_idx[i] is not None: #if nearby ESL data found, try to analyze that data
    
            matched_filenames.append([station_filenames[x] for x in min_idx[i]]) #esl stations within distance to sea-level projection site
            sites_with_esl.append(input_locations.locations[i].values) #sea-level projections sites with nearby ESL data
        else:
            print("Warning: No nearby ESL information found for {0}. Skipping site.".format(input_locations.locations[i].values))
            
    # if no locations with esl data are found, quit with error
    if not matched_filenames:
        raise Exception("Zero matches found within {} km for provided lat/lon list".format(match_dist_limit))
      
    # Initialize variables to track the files that have been tested
    pass_files = {}
    fail_files = []
    esl_statistics = {}	
    
    if preproc_settings['store_esls']:
        if os.path.exists(os.path.join(output_dir,'ESLs')) == False:
            os.mkdir(os.path.join(output_dir,'ESLs'))

    # Loop over the matched files
    print('Analyzing GESLA records.')
    for i in tqdm(np.arange(len(sites_with_esl))): #loop over sea-level projection sites with matched ESL data
        # This ID
        this_id = str(sites_with_esl[i])
        this_id_passed = False

        # Loop over the esl files within the match radius for this location, from nearest to furthest, use the first one for which data fulfills criteria
        
        ### an alternative option could be to open all fns within radius, and use the longest?
        for esl_file in matched_filenames[i]:
            
            # if this esl file was already tested and passed
            if np.isin(esl_file, list(pass_files.keys())):
                print("{0} already PASSED a previous check on data constraints. Mapping site ID {1} to site ID {2}.".format(esl_file, pass_files[esl_file], this_id))
                esl_statistics[this_id] = esl_statistics[pass_files[str(esl_file)]]
                this_id_passed = True
            elif np.isin(esl_file,fail_files):
                continue
    
            # if current sea-level projection location already has esl parameters, skip to the next
            if this_id_passed:
                continue            
            
            # This esl file has not been tested yet, try to analyze it
            try:
                #this tide gauge data
                if gesla_version == 'gesla2':
                    dfs = ingest_GESLA2_files(path_to_gesla,preproc_settings,fns=[esl_file])
    
                elif gesla_version == 'gesla3':
                    dfs = ingest_GESLA3_files(path_to_gesla,preproc_settings,fns=[esl_file])
                    
                ###   
                #preproc options:
                if preproc_settings['detrend']:
                    dfs = detrend_gesla_dfs(dfs)
                if preproc_settings['deseasonalize']:
                    dfs = deseasonalize_gesla_dfs(dfs)
                if preproc_settings['subtract_amean']:
                    dfs = subtract_amean_from_gesla_dfs(dfs)
                
                extremes = pot_extremes_from_gesla_dfs(dfs,preproc_settings['extremes_threshold'],preproc_settings['declus_method'],preproc_settings['declus_window'])
                
                if preproc_settings['store_esls']:
                    for k,v in extremes.items():
                        v.to_csv(os.path.join(output_dir,'ESLs',k+'_extremes.csv'))
                        
                gpd_params = fit_gpd_to_gesla_extremes(extremes)
                
                esl_statistics[this_id] = gpd_params
                esl_statistics[this_id]['vdatum'] = dfs[esl_file].attrs['vdatum']
                esl_statistics[this_id]['matched_lon'] = dfs[esl_file].attrs['longitude']
                esl_statistics[this_id]['matched_lat'] = dfs[esl_file].attrs['latitude']
                esl_statistics[this_id]['matched_id'] = esl_file
                # This file passed, add it to the pass file dictionary
                pass_files[str(esl_file)] = this_id
                # This ID has data, set the flag to move onto the next ID
                this_id_passed = True
    			
                #possibly future work: code an optional storage of the extremes here, probabily useful
                
            except:
                # This file failed, add it to the fail list and continue with the next file
                fail_files.append(str(esl_file))
                continue
    
        # Let the user know we didn't find a file that passes the data constraints
        if not this_id_passed:
            print("No locations within {0} km pass the data constraints for ID {1}.".format(match_dist_limit, this_id))
        
    if len(esl_statistics) == 0:
        raise Exception('Did not find any nearby records passing the data constraints.')
    return esl_statistics

def pot_extremes_from_gesla_dfs(dfs,threshold_pct,declus_method=None,declus_window=None):
    ''' Derive declustered peaks over threshold based on preprocessed data in "dfs". Declustering is optional but recommended.'''
    
    assert ((threshold_pct>0) & (threshold_pct<100))
    extremes_dfs = {}
    
    for k,df in dfs.items():
        threshold = np.percentile(df['sea_level'],threshold_pct)
        extremes = df[df['sea_level']>=threshold]

        
        if declus_method == 'iterative_descending':
            #sorts extremes descendingly, and iteratively add extremes to declustered extremes if they are further than declus_window away from already included extremes
            extremes_sorted = extremes.sort_values(by=['sea_level'], ascending=False, kind='mergesort')
            
            decl_extremes = extremes_sorted.iloc[[0]] #initialize using the first (highest) extreme
            
            for i in range(1,len(extremes_sorted)): #for each subsequent extreme
                next_extreme = extremes_sorted.iloc[[i]]
                timediffs = decl_extremes.index - next_extreme.index[0] #calculate time differences
            
                # if none of the included extremes are less than 'window_len' days away from current extreme, add current extreme to declustered extremes
                if all( abs(timediffs / np.timedelta64(24*3600, 's')) >= declus_window): 
                    decl_extremes = pd.concat((decl_extremes,next_extreme)) #add next extreme
            
            extremes_dfs[k] = decl_extremes.sort_index()

        elif declus_method == 'rolling_max':
            #only include extremes that are the maximum of a rolling window of declus_window around that extreme
            extremes_dfs[k] = extremes['sea_level'].where(extremes['sea_level']==extremes['sea_level'].rolling(window=str(declus_window)+'D',center=True,min_periods=1).max()).dropna().to_frame()
            
        elif declus_method == None:
            extremes_dfs[k] = extremes
        
        if 'H' in df.attrs['resample_freq']:
            rps = pd.DataFrame(data=[(1 + len(df)/365.25/24)/ i for i in range(1,len(extremes_dfs[k])+1)],
                               index=extremes_dfs[k].sort_values(by=['sea_level'], ascending=False, kind='mergesort').index,
                               columns=['rp_empirical']) #determine return periods of declustered extremes empirically
            extremes_dfs[k] = extremes_dfs[k].join(rps.sort_index())
            
            extremes_dfs[k].attrs['avg_extr_pyear'] = 24 * 365.25 * len(extremes_dfs[k])/len(df)
        elif 'D' in df.attrs['resample_freq']:
            rps = pd.DataFrame(data=[(1 + len(df)/365.25)/ i for i in range(1,len(extremes_dfs[k])+1)],
                               index=extremes_dfs[k].sort_values(by=['sea_level'], ascending=False, kind='mergesort').index,
                               columns=['rp_empirical'])
            extremes_dfs[k] = extremes_dfs[k].join(rps.sort_index())
            
            extremes_dfs[k].attrs['avg_extr_pyear'] = 365.25 * len(extremes_dfs[k])/len(df)
        
        #add metadata
        extremes_dfs[k].attrs['mhhw'] = np.nanmean(df.groupby(pd.Grouper(freq='2D')).max()) #2d because some places don't see high tide every day
        extremes_dfs[k].attrs['threshold_pct'] = threshold_pct
        extremes_dfs[k].attrs['threshold'] = threshold
        extremes_dfs[k].attrs['declus_method'] = declus_method
        extremes_dfs[k].attrs['declus_window'] = declus_window
        
    return extremes_dfs

def guess_initial_gpd_params(data):
    '''provide initial guess of GPD parameters for extremes input data using method of moments (code based on MATLAB script gpfit.m)'''
    xbar = np.mean(data) #mean of extremes
    s2 = np.var(data) #variance of extremes
    k0 = -.5 * ((xbar**2)/s2 - 1) #initial guesses
    sigma0 = .5 * xbar * ( (xbar**2) / s2 + 1)
    xmax = max(data)
    
    if (k0 < 0 and xmax >= -sigma0/k0 ):# if method of moments invalid (code based on MATLAB script gpfit.m), #assume exponential distribution
        k0 = 0
        sigma0 = xbar
    
    return k0,sigma0

def fit_gpd_to_extremes_ds(extremes_minus_loc):
    extremes_minus_loc = extremes_minus_loc[np.isfinite(extremes_minus_loc)]

    k0,sigma0 = guess_initial_gpd_params(extremes_minus_loc)
    
    gp_params = genpareto.fit(extremes_minus_loc,loc=0,scale=sigma0,floc=0) #fit gpd based on ESLs and initial guess, note that optimization vals differ slightly from using gpfit.m
    gp_nlogl, gp_cov = gplike(gp_params[0], gp_params[2],extremes_minus_loc) #calculate covariance matrix of estimated parameters
    
    return np.array([gp_params[-1],gp_params[0]]),gp_nlogl,gp_cov #scale/shape params, nlogl MLE diagnostic, scale/shape cov

def fit_gpd_to_gesla_extremes(dfs):
    '''
    For each dataframe in dict of dataframes, do MLE of GPD parameters based on peaks over threshold.
    Estimates are initialized using method of moments and estimated with python version of gplike.m. 
    
    Expected inputs are extremes derived from GESLA data stored under ['sea_level'].
    '''
    for k,df in dfs.items():
        loc = df.attrs['threshold'] #location parameter = threshold
        
        k0,sigma0 = guess_initial_gpd_params(df['sea_level'].values-loc)
    
        gp_params = genpareto.fit(df['sea_level'].values-loc,loc=0,scale=sigma0,floc=0) #fit gpd based on ESLs and initial guess, note that optimization vals differ slightly from using gpfit.m #TH edit 8-2-24: may have to do with not fixing floc=0 as in gpfit.m, changed this
        gp_nlogl, gp_cov = gplike(gp_params[0], gp_params[2],df['sea_level'].values-loc) #calculate covariance matrix of estimated parameters
    
        if gp_nlogl == np.Inf:
            #GPD parameters not supported, confidence intervals and standard errors cannot be computed reliably
            print(k+' - not included, GPD parameter estimation is unreliable.')
            continue
        
        d = {'loc':[loc], 'scale': [gp_params[-1]], 'shape': [gp_params[0]], 'cov': [gp_cov], 'avg_extr_pyear': df.attrs['avg_extr_pyear'], 'mhhw': df.attrs['mhhw'], 'key': [k]}
        gpd_df = pd.DataFrame(data=d).set_index('key')

        if k == list(dfs.keys())[0]:
            gpd_dfs = gpd_df    
        else:
            gpd_dfs = pd.concat((gpd_dfs,gpd_df))
            
    return gpd_dfs #output dictionaries with GPD parameters


def gplike(shape,scale,data):
    """ gplike.py
    Python version of MATLAB's 'gplike':
        
        #GPLIKE Negative log-likelihood for the generalized Pareto distribution.
        #NLOGL = GPLIKE(PARAMS,DATA) returns the negative of the log-likelihood for
        #the two-parameter generalized Pareto (GP) distribution, evaluated at
        #parameters shape = K and scale = SIGMA, given DATA. GPLIKE does
        #not allow a threshold (location) parameter. NLOGL is a scalar.
        #
        #[NLOGL, ACOV] = GPLIKE(PARAMS,DATA) returns the inverse of Fisher's
        #information matrix, ACOV.  If the input parameter values in shape, scale are the
        #maximum likelihood estimates, the diagonal elements of ACOV are their
        #asymptotic variances.  ACOV is based on the observed Fisher's information,
        #not the expected information.
        #
        #When K = 0 and THETA = 0, the GP is equivalent to the exponential
        #distribution.  When K > 0 and THETA = SIGMA/K, the GP is equivalent to the
        #Pareto distribution.  The mean of the GP is not finite when K >= 1, and the
        #variance is not finite when K >= 1/2.  When K >= 0, the GP has positive
        #density for X>THETA, or, when K < 0, for 0 <= (X-THETA)/SIGMA <= -1/K.
        #
        #See also GPCDF, GPFIT, GPINV, GPPDF, GPRND, GPSTAT.
        
        #References:
        #      [1] Embrechts, P., C. Kl�ppelberg, and T. Mikosch (1997) Modelling
        #          Extremal Events for Insurance and Finance, Springer.
        #      [2] Kotz, S. and S. Nadarajah (2001) Extreme Value Distributions:
        #          Theory and Applications, World Scientific Publishing Company.
        
    Parameters:
    shape = shape coefficient of GPD
    scale = scale coefficient of GPD
    data  = data to which the GPD is fitted

    Output: nlogl and acov, see above

    Created on Tue Nov  5 12:24:26 2019
    @author: Tim Hermans
    """
    k       = shape   # Tail index parameter
    sigma   = scale   # Scale parameter
    lnsigma = np.log(sigma);   # Scale parameter, logged
    
    n = len(data);
    z = data/sigma;
    
    if abs(k) > np.spacing(1):
        if k > 0 or max(z) < -1/k:
            u = 1 + k*z
            sumlnu = sum(np.log1p(k*z));
            nlogL = n*lnsigma + (1+1/k)*sumlnu;
            v = z/u;
            sumv = sum(v);
            sumvsq = sum(v**2);
            nH11 = 2*sumlnu/k**3 - 2*sumv/k**2 - (1+1/k)*sumvsq;
            nH12 = (-sumv + (k+1)*sumvsq)/sigma;
            nH22 = (-n + 2*(k+1)*sumv - k*(k+1)*sumvsq)/sigma**2;
            acov = [[nH22, -nH12],[-nH12, nH11]] / (nH11*nH22 - nH12*nH12)
        else:
            # The support of the GP when k<0 is 0 < y < abs(sigma/k)
            nlogL = np.Inf;
            acov = [[np.nan, np.nan], [np.nan, np.nan]];
            
    else: # limiting exponential dist'n as k->0
        # Handle limit explicitly to prevent (1/0) * log(1) == Inf*0 == NaN.
        nlogL = n*lnsigma + sum(z);
        sumz = sum(z);
        sumzsq = sum(z**2);
        sumzcb = sum(z**3);
        nH11 = (2/3)*sumzcb - sumzsq;
        nH12 = (-n + 2*sumz)/sigma**2;
        nH22 = (-sumz + sumzsq)/sigma;
        acov = [[nH22, -nH12], [-nH12, nH11]] / (nH11*nH22 - nH12*nH12)
    
    return nlogL,acov

def infer_avg_extr_pyear(loc,scale,shape,rz,rf):
    '''infers average rate of exceedance of location parameter given GPD parameters and known return curve'''
    return rf/np.power((1+(shape*(rz-loc)/scale)),(-1/shape))


def multivariate_normal_gpd_samples_from_covmat(scale,shape,cov,n,seed=None):
    '''generate scale and shape samples assuming a multivariate normal distribution with covariance matrix "cov"'''
    if seed:
        np.random.seed(seed)
    gp_samples = np.random.multivariate_normal([shape,scale], cov,size=n)
    shape_samples = gp_samples[:,0]
    scale_samples = gp_samples[:,1]
    scale_samples[scale_samples<0] = 0.00001 #no negative scales, assign very small scale instead (completely 0 leads to division problems)
    
    return scale_samples, shape_samples

def gum_amax_Z_from_F(scale,loc,f): 
    '''compute return heights for return frequencies "f" based on Gumbel distribution parameters'''
    try: #convert array input to float if using only 1 scale
        scale=scale.item()
    except:
        pass
    f = np.where(f<1,f,np.nan) #equation supported to f=1?
    
    z = 0*f #initialize z
    z = loc - np.log( -np.log(1-f) ) * scale #see Eq 3.4 in Coles; inverted from cdf formula exceedeance prob = 1-np.exp(-np.exp((location-z)/scale)
    
    return z 
    
def gpd_Z_from_F(scale,shape,loc,avg_exceed,f): 
    """
    Obtain Z for F, with Z evaluating to NaN for F lower than the average 
    exceedance of the location parameter. Takes into account the lower and 
    upper supported values of the GPD distribution.
    
    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        f:   sample frequencies to compute Z for
        
    Output:
        z:              heights computed for sample frequencies, relative to location parameter
    """
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
    
    z = 0*f #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(f/avg_exceed) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape==0):
            z[shape==0] = -scale[shape==0] * np.log(f[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (f[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 )
            
    #test lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc))
    z = np.where(f>avg_exceed,np.nan,z) #where F>avg_exceed, replace z by np.nan
    
    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    
    #return z #output relative to location parameter
    return z + loc #output relative to vertical datum
    
def gpd_Z_from_F_mhhw(scale,shape,loc,avg_exceed,f,mhhw): 
    """
    Obtain Z for F, with Z for F lower than the average exceedance of the 
    location parameter evaluated using a Gumbel distribution between
    the location parameter and MHHW (with frequency of once per two days). 
    Takes into account the lower and upper supported values of the GPD distribution.
    
    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        f:   sample frequencies to compute Z for
        mhhw:           z corresponding to MHHW (absolute?)
        mhhw_freq:      frequency corresponding to MHHW
    Output:
        z:              heights computed for sample frequencies, relative to location parameter 
                        (meaning z is negative for sample_freq>avg_exceed, i.e., in Gumbel distribution below location parameter)
    """
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
        
    z = 0*f #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(f/avg_exceed) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape==0):
            z[shape==0] = -scale[shape==0] * np.log(f[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (f[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 )

    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    
    #below lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc)), use Gumbel distribution
    z[f>avg_exceed] = 0 + np.log(f[f>avg_exceed]/avg_exceed) * (mhhw-loc)/np.log((365.25/2)/avg_exceed) #lower bound of gumbel is z=z0-loc=(mhhw-loc), i.e. (loc-mhhw) below loc if mhhw follows from long-term mean of 2-day maximum

    #below lower bound of Gumbel?
    z = np.where(f>(365.25/2),np.nan,z) #where F>mhhwFreq, replace z by np.nan

    #return z #output relative to location parameter
    return z + loc #output relative to vertical datum


def gpd_Z_from_F_Sweet22(scale,shape,loc,avg_exceed,f): 
    """
    Obtain Z for F, with Z for F lower than the average exceedance of the 
    location parameter evaluated through extrapolation between .5/yr and the 
    average exceedance, up to 10/yr, following (Sweet et al., 2022).

    Input: 
        scale:          scale parameter (scalar or vector)
        shape:          shape parameter (scalar or vector)
        loc:            location parameter
        avg_exceed:     exceedance frequeny of location parameter (POT threshold)
        f:   sample frequencies to compute Z for
       
    Output:
        z:              heights computed for sample frequencies, relative to location parameter 
                        (meaning z is negative for sample_freq>avg_exceed, i.e., in Gumbel distribution below location parameter)
    """
    
    try: #convert array input to float if using only 1 shape & scale
        scale=scale.item()
        shape=shape.item()
    except:
        pass
    
    z = 0*f #initialize z
    
    if np.isscalar(shape): #if shape is scalar value
        if shape ==0: 
            z = -scale * np.log(f/avg_exceed) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 ) #rearranged from F = f(Z), as in e.g., Frederikse et al. (2020), Buchanan et al. (2016)
            
    else: #if array input
        if np.any(shape==0):
            z[shape==0] = -scale[shape==0] * np.log(f[shape==0]/avg_exceed)
            z[shape!=0] = scale[shape!=0]/shape[shape!=0] * (np.power( (f[shape!=0]/avg_exceed), (shape[shape!=0]/-1)) -1 )   
        else:
            z = scale/shape * (np.power( (f/avg_exceed), (shape/-1)) -1 )

    #test upper bound of GPD ((shape*z/scale)<-1)
    z = np.where((shape*z/scale)<-1,np.nan,z)
    
    #below lower bound of GPD (F=avg_exceed for z->0 with z=z0-loc (z0->loc)), use extrapolation
    z = np.where(f>avg_exceed,np.nan,z)

    z_0p5 = gpd_Z_from_F(scale,shape,loc,avg_exceed,np.array([0.5]))-loc #heights corresponding to 0.5/yr frequency
    z_0p2 = gpd_Z_from_F(scale,shape,loc,avg_exceed,np.array([0.2]))-loc #heights corresponding to 0.2/yr frequency
    
    iXtrp = ((f>avg_exceed) & (f<=10)) #indices of frequencies to fill with extrapolation
    
    if ((avg_exceed>0.5)&(avg_exceed<10)):
        if np.isscalar(shape): #if shape is scalar value
            z[iXtrp] = z_0p5+np.log(f[iXtrp]/0.5) * (0-z_0p5)/np.log(avg_exceed/0.5) #log-linear similar to Gumbel
        else:
            
            z[iXtrp] = z_0p5[iXtrp]+np.log(f[iXtrp]/0.5) * (0-z_0p5[iXtrp])/np.log(avg_exceed/0.5)
    
    elif ((avg_exceed>0.2)&(avg_exceed<10)):
        if np.isscalar(shape): #if shape is scalar value
            z[iXtrp] = z_0p2+np.log(f[iXtrp]/0.2) * (0-z_0p2)/np.log(avg_exceed/0.2) #log-linear similar to Gumbel
        else:
            z[iXtrp] = z_0p2[iXtrp]+np.log(f[iXtrp]/0.2) * (0-z_0p2[iXtrp])/np.log(avg_exceed/0.2)
    
    if avg_exceed<=0.2:
        print('Warning: "avg_extr_pyear" is equal to or less than 0.2/yr. Cannot extrapolate using the method of Sweet et al. (2022).')
    #return z #output relative to location parameter
    return z + loc #output relative to vertical datum


def get_return_curve_gpd(f,scale,shape,loc,avg_exceed,below_gpd=None,mhhw=None):
    '''wrapper function around gpd_Z_from_F functions that handles input/output shapes'''
    assert np.shape(scale)==np.shape(shape)
    
    f_=f
    if np.isscalar(scale) == False:
        f_ = np.repeat(f_[:,np.newaxis],len(scale),axis=1) #repeat f num_mc times
        scale = np.repeat(scale[np.newaxis,:],len(f),0)
        shape = np.repeat(shape[np.newaxis,:],len(f),0)
        
    z=np.nan*f_ #initialize
    
    #compute z from f using gpd parameters
    if below_gpd == 'mhhw':
        z = gpd_Z_from_F_mhhw(scale,shape,loc,avg_exceed,f_,mhhw)
        
    elif below_gpd == 'Sweet22':
        z = gpd_Z_from_F_Sweet22(scale,shape,loc,avg_exceed,f_)
    
    elif not below_gpd:
        z = gpd_Z_from_F(scale,shape,loc,avg_exceed,f_)
    
    return z

def get_return_curve_gumbel(f,scale,loc):
    '''wrapper function around gum_Z_from_F functions that handles input/output shapes'''
    f_=f
    if np.isscalar(scale) == False:
        f_ = np.repeat(f_[:,np.newaxis],len(scale),axis=1) #repeat f num_mc times
        scale = np.repeat(scale[np.newaxis,:],len(f),0) 
        
    z=np.nan*f_ #initialize
    
    #compute z from f using gpd parameters
    z = gum_amax_Z_from_F(scale,loc,f_)
    
    return z