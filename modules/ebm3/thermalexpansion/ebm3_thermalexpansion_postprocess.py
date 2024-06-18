import numpy as np
import os
import pandas as pd
import sys
import argparse
from fair import energy_balance_model as ebm3
import csv
import xarray as xr
from sklearn import linear_model
import re
import time
import netCDF4
import h5py
import scipy

"""
Created on Fri Jun  7 11:24:19 2024

For any model you have up tp 2300, you want a set of scaling coefficients 
for ssp126 and ssp585, as those thwo a typically the ones available up to 2300.
Those scenarios are used to get a pattern of regression coefficient used for
SSP interpolation and probabilistic projections using FAIR.

@author: vmalagonsantos
"""

# functions 
def fetch_erfs_from_rcmip(path,scenarios): 
    '''  
    function to read scenario-dependent ERF from RCMIP csv
    
    Input parameters: 
    path        = path to RCMIP csv
    scenarios   = desired ssp's

    Output: Adds results of extreme sea-level analysis for all stations to the 
    station data dictionary.
    
    erfs        = effective radiative forcing per scenario and year
    erfyears    = corresponding years (1750-2500)
    '''

    erfyears=np.arange(1750,2501)
    ssp_idx = {'ssp119':212,'ssp126':231,'ssp245':308,'ssp370':59,'ssp585':404} #location in table
    
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file)
        rows = list(csv_reader)
    
    erfs = np.empty((len(scenarios),len(erfyears)))
    
    for s,scen in enumerate(scenarios):
        try:    
            erfs[s,:] = rows[ssp_idx[scen]][7:]
        except:
            continue
        
    return erfs, erfyears

def angd(lat0, lon0, qlat, qlon):

	# Convert the input from degrees to radians
	(lat0, lon0) = np.radians((lat0, lon0))
	(qlat, qlon) = np.radians((qlat, qlon))

	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(qlat)*np.sin(qlon-lon0))**2 + \
	(np.cos(lat0)*np.sin(qlat) - np.sin(lat0)*np.cos(qlat) * np.cos(qlon-lon0))**2),\
	(np.sin(lat0)*np.sin(qlat) + np.cos(lat0)*np.cos(qlat)*np.cos(qlon-lon0)))

	# Convert the results from radians to degrees and return
	return(np.degrees(temp))

def NearestPoint(qlat, qlon, lats, lons, tol = None):

	# Get the distance between the query point and all the possible points
	dist = angd(lats, lons, qlat, qlon)

	# Which is the closest point
	nearest_idx = np.argmin(dist)

	# Is the point within the tolerance?
	if isinstance(tol, (int, float)):
		if dist[nearest_idx] > tol:
			return(None)

	return(nearest_idx)


def NearestPoints(qlats, qlons, lats, lons, tol):

	if len(qlats) != len(qlons):
		raise Exception("Query lats ({}) and lons ({}) differ in length".format(len(qlats), len(qlons)))

	idx = map(lambda qlat,qlon: NearestPoint(qlat, qlon, lats, lons, tol), qlats, qlons)

	return(list(idx))

def ReadLocationFile(location_file):

	# Initialize variables to hold data and site information
	names = []
	ids = []
	lats = []
	lons = []
	
	# Compile the regex for finding commented lines
	comment_regex = re.compile(r'^#')
	
	# Open the rate file
	with open(location_file, 'r') as f:
		
		# Loop over the lines of the file
		for line in f:
			
			# Skip commented lines
			if re.search(comment_regex, line):
				continue
			
			# Split the line into components
			(this_name, this_id, this_lat, this_lon) = line.split("\t")

			# Store the information
			names.append(this_name)
			ids.append(int(this_id))
			lats.append(float(this_lat))
			lons.append(float(this_lon))
	
	# Cast everything as numpy arrays
	names = np.array(names)
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	
	# Return variables
	return(names, ids, lats, lons)


# Smooth ZOS and ZOSTOGA over 19 year smoothing window
def nanSmooth(x, w=19):
	idx = np.flatnonzero(~np.isnan(x))
	temp = x
	if len(idx) > 0:
		temp[idx] = Smooth(x[idx], w)
	return(temp)

def Smooth(x, w=19):
	out0 = np.convolve(x, np.ones(w,dtype='double'), 'valid')/w
	r = np.arange(1,w-1,2, dtype="double")
	start = np.cumsum(x[:w-1])[::2]/r
	stop = (np.cumsum(x[:-w:-1])[::2]/r)[::-1]
	y = np.concatenate((start, out0, stop))
	return(y)

def emb3_thermalexpansion_postprocess(scenario, pipeline_id, nsamps, seed, pyear_start, pyear_end, pyear_step, locationfile, baseyear, climate_data_file, rfmip, params):
    targyears = np.arange(baseyear,2301) # for regression
    projyears = np.arange(pyear_start,pyear_end+1, pyear_step)
    (_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

    # get temperature from FaIR simulations: INPUT FROM CLIMATE STEP. We need both gsat and oceantemp
    cf = xr.open_dataset(climate_data_file, group=scenario, engine='netcdf4', )
    gsat = cf['surface_temperature'].sel(years=projyears) - cf['surface_temperature'].sel(years=np.arange(baseyear-9,baseyear+10)).mean(dim='years')
    otemp = cf['deep_ocean_temperature'].sel(years=projyears) - cf['deep_ocean_temperature'].sel(years=np.arange(baseyear-9,baseyear+10)).mean(dim='years')

    # INPUT temperature file from ebm3 global
    gte_file = 'ssp585.emu2.2300.fair2.ocean.tlm.sterodynamics_globalsl.nc'

    # INPUT get models and parameters
    zosdir = 'cmip6/zos/' # INPU cmip6 zos
    paramdir = params

    #forcing for ebm
    scenarios = ['ssp126', 'ssp585'] # NOT INPUT, these two are needed for SSP interpolation in 2300 projections
    rfmipfile = rfmip #path # RFMIP FILE IS INPUT, NEW DATA
    erfs,erfyears = fetch_erfs_from_rcmip(rfmipfile, scenarios) #get ERF timeseries

    ebm_param = pd.read_csv(paramdir)

    # find available zos simulations - ignoring hidden files
    model_has_zos = [f for f in os.listdir(zosdir) if not f.startswith('.')]

    ## FITTING
    slopes = []
    intercept = []

    for m, model in enumerate(model_has_zos):
    # model = model_has_zos[0]

        # model = 'MRI-ESM2-0'
        
        # find matching ebm parameters to model
        model_param = ebm_param[ebm_param["model"].str.contains(model)]
        
        modeldir = zosdir + model
        # keep those up to 2300 or 2500
        zos_runs2_2300 = [f for f in os.listdir(modeldir) if '2300' in f or '2500' in f] # regex would be better here
        if not zos_runs2_2300:
            print('Model ' + model + ' only runs to 2100. Skipping and looking for longer simulations' )

            continue
        
        print('Estimating parameters for model: ' + model)

            
        # get variants with 2300 simulation
        if len(zos_runs2_2300) == 1:
            variant_2300 = [f for f in re.split('_', zos_runs2_2300[0]) if f.startswith('r')]
            model_var_param = model_param[ebm_param["run"].str.contains(variant_2300[0])]
        else:
            for v, filename in enumerate(zos_runs2_2300):
                variant_2300 = [f for f in re.split('_', zos_runs2_2300[v]) if f.startswith('r')]
                # get EBM parameters for both model and variant
                model_var_param = model_param[ebm_param["run"].str.contains(variant_2300[0])]
                if not model_var_param.empty:
                    continue
        
        c1 = np.array(model_var_param.get('C1'))
        c2 = np.array(model_var_param.get('C2'))
        c3 = np.array(model_var_param.get('C3'))
        k1 = np.array(model_var_param.get('kappa1'))
        k2 = np.array(model_var_param.get('kappa2'))
        k3 = np.array(model_var_param.get('kappa3'))
        e = np.array(model_var_param.get('epsilon'))
        f4 = np.array(model_var_param.get('F_4xCO2'))
        se = np.array(model_var_param.get('sigma_eta'))
        sx = np.array(model_var_param.get('sigma_xi'))
        gamma = np.array(model_var_param.get('gamma'))
        nit = np.array(model_var_param.get('nit'))
        
        
        temp3 = np.empty((len(erfyears), 3, len(scenarios)))
        
        for s, scen in enumerate(scenarios):
        
        #Model setup 3lm
            ebm_obj = ebm3.EnergyBalanceModel(
                    ocean_heat_capacity = [c1[0], c2[0], c3[0]],
                    ocean_heat_transfer =[k1[0], k2[0], k3[0]], ## change this placeholder value
                    deep_ocean_efficacy = e[0],
                    forcing_4co2 = f4[0],
                    stochastic_run=False,
                    sigma_eta = se[0],
                    sigma_xi = sx[0],
                    gamma_autocorrelation = gamma[0],
                    seed = None,
                    timestep = 1,
                    n_timesteps = 2500-1750,
                )
            
        # Run model for different scenarios
            ebm_obj.add_forcing( erfs[s,:], erfyears)
            ebm_obj.run()
            temp3[:,:,s] = ebm_obj.temperature
        
        attrs={'description': 'surface and deep temperature obtained by EBM-3LM',
                    'contact':'Victor Malagon Santos, victor.malagon.santos@nioz.nl'}
        
        temp3_xr = xr.Dataset({'temperature':(['years','layers', 'scenarios'], temp3)},
        coords={'layers':['surface', 'inter', 'deep'], 'years': erfyears, 'scenarios': scenarios}, attrs=attrs)
        
        # reference to base year
        temp3_xr = temp3_xr - temp3_xr.sel(years=baseyear)
        # .sel(years=np.arange(baseyear-9,baseyear+10)).mean(dim='years')
        
        # get zos, historical + ssp
        model_dir = zosdir + model
        
        hist_file = [f for f in os.listdir(model_dir) if 'hist' in f and variant_2300[0] in f]
        ssp1_file = [f for f in os.listdir(model_dir) if 'ssp1' in f and variant_2300[0] in f and 'ssp126' in f]
        ssp5_file = [f for f in os.listdir(model_dir) if 'ssp5' in f and variant_2300[0] in f and 'ssp585' in f]
        
        if not ssp1_file or not ssp5_file:
            continue
        
        hist_xr = xr.open_dataset(model_dir + '/' + hist_file[0],decode_times=False)
        ssp1_xr = xr.open_dataset(model_dir + '/' + ssp1_file[0],decode_times=False)
        ssp5_xr = xr.open_dataset(model_dir + '/' + ssp5_file[0],decode_times=False)
        
        if ssp1_xr['time'][0].values != 60280.5:
            ssp1_xr.coords['time'] = ((ssp1_xr['time'] + 60280.5).astype('float'))
            ssp5_xr.coords['time'] = ((ssp5_xr['time'] + 60280.5).astype('float'))
        
        # cocatenate
        zos1 = xr.concat([hist_xr, ssp1_xr], dim='time')
        zos5 = xr.concat([hist_xr, ssp5_xr], dim='time')
        
        #% create new coordinates for year
        zos1.coords['years'] = ((zos1['time'] // 365.25)+1850).astype('int')   # REMOVE HARDCODED STARTING YEAR
        zos5.coords['years'] = ((zos5['time'] // 365.25)+1850).astype('int')   # REMOVE HARDCODED STARTING YEAR
        
        # Calculate annual means
        zos1_am = zos1.groupby('years').mean(dim='time')
        zos1_am = zos1_am.where(zos1_am['zos'] < 99999) # replace land values with nans
        zos5_am = zos5.groupby('years').mean(dim='time')
        zos5_am = zos5_am.where(zos5_am['zos'] < 99999) # replace land values with nans

        # reference zos
        # baseyear_idx = np.flatnonzero(datayears == baseyear)
        # S = np.apply_along_axis(lambda z, idx: z - z[idx], axis=0, arr=sZOS, idx=baseyear_idx)
        zos1_am = zos1_am.sel(years=targyears) - zos1_am.sel(years=baseyear)
        zos5_am = zos5_am.sel(years=targyears) - zos5_am.sel(years=baseyear)
        
        #% fitting stage
        nlon = len(zos1_am['lon'])
        nlat = len(zos1_am['lat'])
            
        
        Ts = xr.concat((temp3_xr.sel(years=targyears).sel(scenarios=['ssp126']).sel(layers='surface').drop_vars('scenarios'),
                        temp3_xr.sel(years=targyears).sel(scenarios=['ssp585']).sel(layers='surface').drop_vars('scenarios')),dim='years')
        
        Ti = xr.concat((temp3_xr.sel(years=targyears).sel(scenarios=['ssp126']).sel(layers='inter').drop_vars('scenarios'),
                        temp3_xr.sel(years=targyears).sel(scenarios=['ssp585']).sel(layers='inter').drop_vars('scenarios')),dim='years')
        
        Td = xr.concat((temp3_xr.sel(years=targyears).sel(scenarios=['ssp126']).sel(layers='deep').drop_vars('scenarios'),
                        temp3_xr.sel(years=targyears).sel(scenarios=['ssp585']).sel(layers='deep').drop_vars('scenarios')),dim='years')
        
        
        zos1 = np.array(zos1_am['zos'].sel(years=targyears))
        zos5 = np.array(zos5_am['zos'].sel(years=targyears))
        
        slope = np.zeros((3, nlat, nlon))
        slope[:] = np.nan
        smoothwin = 19;
        
        reg = linear_model.LinearRegression()
        # for i in tqdm(range(nlat)):
        #     sleep(3)
        for i in tqdm(range(nlat)):
            for j in range(nlon):
        
                y =  np.concatenate((Smooth(zos1[:,i,j],w=smoothwin), Smooth(zos5[:,i,j],w=smoothwin))) # reducing varability
                x1 = np.array(Ts['temperature'][:]).flatten()
                x2 = np.array(Ti['temperature'][:]).flatten()
                x3 = np.array(Td['temperature'][:]).flatten()
                x = np.transpose(np.stack((x1,x2,x3)))
                
                try:
            
                    reg.fit(x, y)
                    slope[0,i,j] = reg.coef_[0]  # slope coefficient of T surface
                    slope[1,i,j] = reg.coef_[1]  # slope coefficient of T inter
                    slope[2,i,j] = reg.coef_[2]  # slope coefficient of T deep
                    # intercept[m,i,j] = reg.intercept_
        
                except ValueError:     
                    slope[0,i,j] = np.nan
                    slope[1,i,j] = np.nan
                    slope[2,i,j] = np.nan
                    # intercept[m,i,j] = np.nan
                    
        slopes.append(slope) # MAYBE SAVE TO A PICKLE, TO BE PROVIDED IN THE NEXT TASK?

    lons, lats = np.meshgrid(zos1_am['lon'],  zos1_am['lat'])
    lats = lats.flatten()
    lons = lons.flatten()

    site_ids_map = np.array(NearestPoints(site_lats, site_lons, lats, lons, tol=None))

    nsamps = 100
    samples = np.array(gsat['samples'])[0:nsamps] # remove preselections of samples

    # RESAMPLE SLOPE PARAMETERS
    nsims=len(slopes)
    nsamps = len(samples)

    rng = np.random.default_rng(seed)
    if nsamps > nsims:
        run_idx = np.arange(nsims)
        sample_idx = rng.choice(nsims, nsamps, nsamps>nsims)
    else:
        run_idx = rng.choice(nsims, nsamps, nsamps>nsims)
        sample_idx = np.arange(nsamps)
        
        
    ## need to figure out of to do this resamples properly, so from now on this is dummy data
    slopes_resampled = [slopes[i] for i in sample_idx]
    # aa = np.dstack(slopes_resampled)

    dsl = np.empty((nsamps, len(projyears), len(site_ids)))
    dsl[:] = np.nan

    for sample in tqdm(samples):
        slope = slopes_resampled[sample]
        print(sample)
        
        # prepare slopes
        slope_s = slope[0,:,:].flatten()
        slope_s = np.array([slope_s[x] for x in site_ids_map])
        
        slope_i = slope[1,:,:].flatten()
        slope_i = np.array([slope_i[x] for x in site_ids_map])
        
        slope_d = slope[2,:,:].flatten()
        slope_d = np.array([slope_d[x] for x in site_ids_map])
        
        # obtain fair temepratures
        Tfs = np.array(gsat['surface_temperature'].sel(locations=-1).sel(samples=sample)) # surface temeprature from fair
        Tfi = np.array(otemp['deep_ocean_temperature'].sel(locations=-1).sel(layers=1).sel(samples=sample)) # intermediate temeprature from fair
        Tfd = np.array(otemp['deep_ocean_temperature'].sel(locations=-1).sel(layers=2).sel(samples=sample)) # deep temeprature from fair
        
        
        #% project
        dsl[sample,:,:] = np.multiply.outer(Tfs, slope_s) + np.multiply.outer(Tfi, slope_i) + np.multiply.outer(Tfd, slope_d) 
        
        # dsl = np.multiply.outer(Tfs, slope_s) + np.multiply.outer(Tfi, slope_i) + np.multiply.outer(Tfd, slope_d) 

    ncvar_attributes = {"description": "Dynamic Sea Level"}

    nc_missing_value = np.nan
    # Generate the output xarray
    dsl_xr = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), dsl, {"units":"mm", "missing_value":nc_missing_value}),
                            "lat": (("locations"), site_lats),
                            "lon": (("locations"), site_lons)},
                            coords={"years": projyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)


    # make sure projections have 0 mean. This weighted average only works for regular grids. 
    # Must be edited to accomodate other grids if zos is not in 1x1.

    weights = np.cos(np.deg2rad(dsl_xr.lat))
    weights.name = "weights"

    dsl_xr_weighted = dsl_xr.weighted(weights)
    dsl_xr_weighted
    weighted_mean = dsl_xr_weighted.mean(("locations"))


    # Appply 0mean correction to dsl, and add GTE to get sterodynamics

    gte = xr.open_dataset(gte_file)


    gte = gte['sea_level_change'].values[:,:,0]
    sdsl = dsl*1000 - np.repeat(np.array(weighted_mean['sea_level_change'])[:, :, np.newaxis], len(site_ids), axis=2) + np.repeat(gte[0:nsamps, :, np.newaxis], len(site_ids), axis=2) 

    ncvar_attributes = {"description": "Local SLR contributions from thermal expansion and dynamic sea-level using EBM3",
            "history": "Created " + time.ctime(time.time()),
            "source": "SLR Framework: PROTECT 2300",
            "scenario": scenario, # change to scenario later
            "baseyear": baseyear}

    nc_missing_value = np.nan
    # Generate the output xarray
    local_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), sdsl, {"units":"mm", "missing_value":nc_missing_value}),
                            "lat": (("locations"), site_lats),
                            "lon": (("locations"), site_lons)},
                            coords={"years": projyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)
        # Write these samples to a temporary netcdf file
    local_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

if __name__ == '__main__':

    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the pre-processing stage for the TLM ocean dynamics workflow",\
    epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

    # Define the command line arguments to be expected
    parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)", default='ssp585')
    parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
    parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
    parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
    parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
    parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=5]", default=5, type=int)
    parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
    parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
    parser.add_argument('--climate_data_file',type=str)
    parser.add_argument('--rfmip', help='rfmip file',default='rfmip-radiative-forcing-annual-means-v4-0-0.csv')
    parser.add_argument('--params', help='CMIP6 Params cvs', default='4xCO2_impulse_response_ebm3_cmip6.csv')

    # Parse the arguments
    args = parser.parse_args()

    emb3_thermalexpansion_postprocess(args.scenario, 
                                      args.pipeline_id, 
                                      args.nsamps, 
                                      args.seed, 
                                      args.pyear_start, 
                                      args.pyear_end, 
                                      args.pyear_step, 
                                      args.locationfile, 
                                      args.baseyear, 
                                      args.climate_data_file,
                                      args.rfmip,
                                      args.params)

    # Done
    sys.exit()