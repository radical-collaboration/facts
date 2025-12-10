""" .py
@author: Tim Hermans
t(dot)h(dot)j(dot)hermans@uu(dot)nl
"""
import numpy as np
import pandas as pd
from gesla import GeslaDataset
import os
from copy import deepcopy
import xarray as xr

def resample_data(rawdata,resample_freq):
    '''Resample raw tide gauge data "rawdata" (dict) to new "resample_freq" (str: ['H_mean','D_mean','D_max'])'''
    h_means = rawdata.resample('H').mean() #as most records provide hourly obs, this is often the mean of 1 obs = the obs
    h_means.loc[h_means['use_flag']!=1,'sea_level'] = np.nan #set hours during which 1 or more observations are bad to nan
    h_means = h_means[np.isfinite(h_means['sea_level'])] #drop nans
    
    if resample_freq == 'H_mean':
        output = h_means
    elif resample_freq == 'D_mean':
        hours_per_day = h_means.resample('D').count()
        d_means = h_means.resample('D').mean()
        d_means = d_means[hours_per_day['sea_level']>=12] #only use days with 12 hourly means or more
        d_means = d_means[np.isfinite(d_means['sea_level'])] #drop nans
        output = d_means
    elif resample_freq == 'D_max':
        hours_per_day = h_means.resample('D').count()
        d_maxs = h_means.resample('D').max()
        d_maxs = d_maxs[hours_per_day['sea_level']>=12] #only use days with 12 hourly means or more
        d_maxs = d_maxs[np.isfinite(d_maxs['sea_level'])] #drop nans
        output = d_maxs
    else:
        raise Exception('Unknown resample frequency.')
        
    output.attrs['resample_freq'] = resample_freq
    return output

def readmeta_GESLA2(filename):
    ''' open GESLA2 header containing metadata'''
    with open(filename,encoding = 'raw_unicode_escape') as myfile:
        head = [next(myfile) for x in range(9)]
    station_name = head[1][12:-1].strip()	 
    station_lat = float(head[4][11:-1].strip())
    station_lon = float(head[5][12:-1].strip()) 
    start_date = np.datetime64(head[7][18:-1].strip().replace('/','-'))
    end_date = np.datetime64(head[8][16:-1].strip().replace('/','-'))
	
    return (station_name, station_lat, station_lon, start_date, end_date)

def extract_GESLA2_locations(path_to_gesla2,min_yrs):
    '''Generate a list of the coordinates of all gesla2 files in "gesladir" '''
    gesladir = path_to_gesla2
    # Get a list of the gesla database files
    geslafiles = os.listdir(gesladir)
    if '.DS_Store' in geslafiles:
        geslafiles.remove('.DS_Store')
    # Initialize the station variables
    station_names = []
    station_lats = []
    station_lons = []
    station_filenames = []
    
    # Loop over the gesla files
    for this_file in geslafiles:
        #Extract the station header information
        this_name, this_lat, this_lon, start_date, end_date = readmeta_GESLA2(os.path.join(gesladir, this_file))
        
        if ((pd.to_datetime(end_date).year - pd.to_datetime(start_date).year) < min_yrs-1) or (pd.to_datetime(end_date).year<2000):
            continue
        # append this information to appropriate lists
        station_names.append(this_name)	 
        station_lats.append(float(this_lat))
        station_lons.append(float(this_lon))
        station_filenames.append(this_file)
    
    return (station_names,station_lats,station_lons,station_filenames)

def extract_GESLA3_locations(path_to_gesla3,min_yrs):
    '''Generate a list of the coordinates of all gesla3 files in using the metadata at "metafile_path" '''
    meta = pd.read_csv(os.path.join(path_to_gesla3,'GESLA3_ALL.csv'))
    
    meta = meta[(meta['GAUGE TYPE']=='Coastal') & (meta['OVERALL RECORD QUALITY']=='No obvious issues') & (meta['NUMBER OF YEARS'] >= min_yrs) & (pd.to_datetime(meta['END DATE/TIME'],format='mixed').dt.year>=2000)].reset_index(drop=True) #only consider coastal gauges without issues with sufficient length
    
    return (list(meta['SITE NAME']),list(meta['LATITUDE']),list(meta['LONGITUDE']),list(meta['FILE NAME']))
    
def get_data_starting_index(file_name):
    '''functionality to find first row with data to read in in GESLA files'''
    with open(file_name,errors='ignore') as fp:
        for i, line in enumerate(fp.readlines()):
                if line.startswith("#"):
                    continue
                else:
                    break
        return i

def reference_data_to_msl(period,data,fn):
    ''' subtract MSL during "MSL_period" from data if possible'''
    y0 = period[0]
    y1 = str(int(period[1])+1)
    
    data_in_msl_period = data[data.index.to_series().between(y0,y1)]
  
    if len(np.unique(data_in_msl_period.index.date)) < 0.5 * 365.25 * (int(y1)-int(y0)):
        print('Warning: Could not reference observations to MSL for {0} because observations are avaible on less than half of all days during "msl_period".'.format(fn))
        vdatum = 'original'
    else:    
        msl = np.nanmean(data_in_msl_period['sea_level'])
        data['sea_level'] = data['sea_level'] - msl
        vdatum = 'MSL (' +y0+'-'+y1+')'
   
    return data, vdatum

def ingest_GESLA3_files(gesla3_path,preproc_settings,fns=None):
    '''open & preprocess files in list "fns" with type in "types" e.g., ['Coastal','River'] (as defined by GESLA) if fulfilling inclusion criteria set in cfg'''
    path_to_files = os.path.join(gesla3_path,'GESLA3.0_ALL')
    path_to_files = os.path.join(path_to_files,'') #append '/'
    
    meta_fn = os.path.join(gesla3_path,'GESLA3_ALL.csv')
    
    resample_freq = preproc_settings['resample_freq']
    min_yrs = preproc_settings['min_yrs']
    
    if not fns:
        fns = os.listdir(path_to_files) #if fns is undefined, try to read in all files
    g3object = GeslaDataset(meta_file=meta_fn,data_path=path_to_files) #create dataset class
    
    datasets = {} #initialize dictionary
    for fn in fns:
        data = g3object.file_to_pandas(fn) #data [0] + metadata [1]
        
        if data[1]['gauge_type'] != 'Coastal':
            continue
        if data[1]['number_of_years'] < min_yrs:
            continue
        
        rawdata = data[0]
    
        if rawdata.index[-1].year < 2000: #if record is older than 2000, skip it
            continue
      
        rawdata.loc[rawdata['use_flag']!=1,'sea_level'] = np.nan
  
        #compute MSL from raw data if referencing to MSL
        if preproc_settings['ref_to_msl']:
            rawdata,vdatum = reference_data_to_msl(preproc_settings['msl_period'].split(','),rawdata,fn)
        else:
            vdatum = 'original'  
        resampled_data  = resample_data(rawdata,resample_freq)
        
        #do not include data if shorter than minimum years of observations
        if (('H' in resample_freq) & (len(resampled_data)/365.25/24 < min_yrs)):
            continue
        if (('D' in resample_freq) & (len(resampled_data)/365.25 < min_yrs)):
            continue
        
        for k,v in data[1].items():
            resampled_data.attrs[k] = v
        resampled_data.attrs['vdatum'] = vdatum

        datasets[fn] = resampled_data.drop(columns=['qc_flag','use_flag'])

    if len(datasets)==0:
        raise Exception('0 records fulfill inclusion criteria.')    
    return datasets


def ingest_GESLA2_files(gesla2_path,preproc_settings,fns=None):
    '''open & preprocess files in list "fns" if fulfilling inclusion criteria set in cfg'''
    path_to_files = gesla2_path
    path_to_files = os.path.join(path_to_files,'') #append '/'
    
    resample_freq = preproc_settings['resample_freq']
    min_yrs = preproc_settings['min_yrs']
    
    if not fns:
        fns = os.listdir(path_to_files) #if fns is undefined, try to read in all files
    
    datasets = {}
    for fn in fns:
        metadata = readmeta_GESLA2(os.path.join(path_to_files,fn))
 
        if (metadata[-1] - metadata[-2])/np.timedelta64(1, 's')/(365.25*24*3600) < min_yrs:
            continue
    
        i = get_data_starting_index(os.path.join(path_to_files,fn))
        
        if 'national_tidal_centre' in fn:
            data = pd.read_csv(
                os.path.join(path_to_files,fn),
                skiprows=i,
                names=['date','time','sea_level','qc_flag','tu','use_flag'],
                sep="\s+",
                parse_dates=[[0, 1]],
                index_col=0,
            )
            
            data = data.drop(columns=['tu'])
        else:
            data = pd.read_csv(
                os.path.join(path_to_files,fn),
                skiprows=i,
                names=['date','time','sea_level','qc_flag','use_flag'],
                sep="\s+",
                parse_dates=[[0, 1]],
                index_col=0,
            )
        if data.index[-1].year < 2000: #if record is older than 2000, skip it
            continue
        
        #set invalid or missing data to nan (see also GESLA2 definitions)
        data.loc[data['use_flag']!=1,'sea_level'] = np.nan
        data.loc[data['sea_level']<=-9.9,'sea_level'] = np.nan
        
        #compute MSL from raw data if referencing to MSL
        if preproc_settings['ref_to_msl']:
            data,vdatum = reference_data_to_msl(preproc_settings['msl_period'].split(','),data,fn)
        else:
            vdatum = 'original'
        resampled_data = resample_data(data,resample_freq)
        
        #do not include data if shorter than minimum years of observations
        if (('H' in resample_freq) & (len(resampled_data)/365.25/24 < min_yrs)):
            continue
        if (('D' in resample_freq) & (len(resampled_data)/365.25 < min_yrs)):
            continue
        
        resampled_data.attrs['file_name'] = fn
        resampled_data.attrs['site_name'] = metadata[0]
        resampled_data.attrs['latitude'] = metadata[1]
        resampled_data.attrs['longitude'] = metadata[2]
        resampled_data.attrs['vdatum'] = vdatum

        datasets[fn] = resampled_data.drop(columns=['qc_flag','use_flag'])
    
    if len(datasets)==0:
        raise Exception('0 records fulfill inclusion criteria.')    
    return datasets

def detrend_gesla_dfs(dfs):
    '''detrend sea_level in dataframes in "dfs"'''
    detrended_dfs = deepcopy(dfs)
    for k,df in detrended_dfs.items():
        x =  df.index.values.astype(np.int64) // 10 ** 9 #convert to seconds timestamp
        y =  df['sea_level'].values.astype('float64')
        lrcoefs = np.polyfit(x,y,1)
        trend = np.polyval(lrcoefs,x)

        df['sea_level'] = df['sea_level'] - trend + np.mean(trend) #subtract trend without changing the vertical datum (mean of timeseries)
    return detrended_dfs

def deseasonalize_gesla_dfs(dfs):
    '''subtract long-term monthly means from sea_level in dataframes in "dfs"'''
    deseasonalized_dfs = deepcopy(dfs)
    for k,df in deseasonalized_dfs.items():

        monthly_means_at_timesteps = df.groupby(df.index.month).transform('mean')['sea_level'].astype('float64') #contains mean of all timesteps in month for all years together at each timestep in that month
        df['sea_level'] = df['sea_level'] - monthly_means_at_timesteps + np.mean(monthly_means_at_timesteps) #subtract without changing the vertical datum (mean of timeseries)
    return deseasonalized_dfs

def detrend_ds(ds,variable):
    p = ds[variable].polyfit(dim='time', deg=1)
    fit = xr.polyval(ds.time, p.polyfit_coefficients)
    
    detrended = ds[variable] - fit
    detrended = detrended + ( ds[variable].mean(dim='time') - detrended.mean(dim='time') )
    
    ds[variable] = detrended
    return ds

def deseasonalize_ds(ds,variable):
    '''subtract long-term monthly means from variable in dataset'''
    deseasoned_ds = ds[variable].groupby(ds.time.dt.month) - ds[variable].groupby(ds.time.dt.month).mean('time')
    
    deseasoned_ds = deseasoned_ds + (ds[variable].mean(dim='time') - deseasoned_ds.mean(dim='time'))
    
    ds[variable] = deseasoned_ds
    return ds

def subtract_amean_from_ds(ds,variable):
    '''subtract annual means from variable in dataset'''
    ds_no_amean = ds[variable].groupby(ds.time.dt.year) - ds[variable].groupby(ds.time.dt.year).mean('time')
    
    ds_no_amean = ds_no_amean + (ds[variable].mean(dim='time') - ds_no_amean.mean(dim='time'))
    
    ds[variable] = ds_no_amean
    return ds

def subtract_amean_from_gesla_dfs(dfs):
    '''subtract annual means from sea_level in dataframes in "dfs"'''
    dfs_no_amean = deepcopy(dfs)
    for k,df in dfs_no_amean.items():
        annual_means_at_timesteps = df.groupby(df.index.year).transform('mean')['sea_level'].astype('float64')
        df['sea_level'] = df['sea_level'] - annual_means_at_timesteps + np.mean(annual_means_at_timesteps)  #subtract without changing the vertical datum (mean of timeseries)
    return dfs_no_amean
    
def drop_shorter_gesla_neighbors(dfs,min_dist=3):
    '''not used in projectESL currently'''
    filtered_dfs = deepcopy(dfs)
    
    for k,df in dfs.items():
        #all tgs in filtered_dfs
        lons    = np.array([k.attrs['longitude'] for k in filtered_dfs.values()])
        lats    = np.array([k.attrs['latitude'] for k in filtered_dfs.values()])
        lengths = np.array([len(k) for k in filtered_dfs.values()])
        
        #current tg
        lon = df.attrs['longitude']
        lat = df.attrs['latitude']
        length = len(df)
        
        #compute distances current tg to all tgs in filtered_dfs
        distances = 6371*2*np.arcsin( np.sqrt(
                np.sin( (np.pi/180) * 0.5*(lats-lat) )**2 +
                np.cos((np.pi/180)*lat)*np.cos((np.pi/180)*lats)*np.sin((np.pi/180)*0.5*(lons-lon))**2) ) #distances from current site to all included sites
        
        if np.sum(distances < min_dist)>1: #if other tgs than current tg within 3 km:
            if length != np.max(lengths[distances<min_dist]): #if current record shorther than nearby records
                filtered_dfs.pop(k) #remove current tg from filtered_dfs
                
    return filtered_dfs