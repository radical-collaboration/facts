'''
@author: Tim Hermans
t(dot)h(dot)j(dot)hermans@uu(dot)nl
'''
import numpy as np
from tqdm import tqdm
from I_O import load_config, open_input_locations, open_gtsm_waterlevels
from preprocessing import extract_GESLA2_locations, extract_GESLA3_locations,ingest_GESLA2_files,ingest_GESLA3_files, detrend_gesla_dfs, deseasonalize_gesla_dfs, subtract_amean_from_gesla_dfs
from preprocessing import detrend_ds, deseasonalize_ds, subtract_amean_from_ds
from esl_analysis import pot_extremes_from_gesla_dfs
from utils import mindist
import os

def save_preprocessed_gesla(gesla_version,path_to_gesla,input_locations,preproc_settings,match_dist_limit,output_dir=None):
    
    if gesla_version == 'gesla2':
        (station_names, station_lats, station_lons, station_filenames) = extract_GESLA2_locations(path_to_gesla,preproc_settings['min_yrs'])
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
        raise Exception("Zero matches found within {} degrees for provided lat/lon list".format(match_dist_limit))
      
    # Initialize variables to track the files that have been tested
    pass_files = {}
    fail_files = []

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
                      
                #preproc options:
                if preproc_settings['detrend']:
                    dfs = detrend_gesla_dfs(dfs)
                if preproc_settings['deseasonalize']:
                    dfs = deseasonalize_gesla_dfs(dfs)
                if preproc_settings['subtract_amean']:
                    dfs = subtract_amean_from_gesla_dfs(dfs)
                
                #create boolean vectors indicating which data exceeds 95th percentile
                extremes = pot_extremes_from_gesla_dfs(dfs,95)
                extremes_declustered = pot_extremes_from_gesla_dfs(dfs,95,'rolling_max',3)
                
                dfs[esl_file]['is_95pct_exceedance'] = 0
                dfs[esl_file]['is_declustered_95pct_exceedance'] = 0
                
                dfs[esl_file].loc[extremes[esl_file].index,'is_95pct_exceedance'] = 1
                dfs[esl_file].loc[extremes_declustered[esl_file].index,'is_declustered_95pct_exceedance'] = 1
                
                dfs[esl_file].to_csv(os.path.join(output_dir,cfg['preprocessing']['resample_freq'] + '_'+esl_file+'.csv'))
                
                pass_files[str(esl_file)] = this_id
                # This ID has data, set the flag to move onto the next ID
                this_id_passed = True
    			
            except:
                # This file failed, add it to the fail list and continue with the next file
                fail_files.append(str(esl_file))
                continue
    return pass_files

def save_preprocessed_gtsm(gtsm_path,input_locations,preproc_settings,match_dist_limit,output_dir=None):
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
        
    threshold= gtsm.waterlevel.quantile(.95,dim='time')
    
    extremes = gtsm.waterlevel.where(gtsm.waterlevel>=threshold)
    gtsm['is_95pct_exceedance'] = np.isfinite(extremes)
    gtsm['is_declustered_95pct_exceedance'] = (extremes==extremes.rolling(time=preproc_settings['declus_window'],center=True,min_periods=1).max()) & (np.isfinite(extremes))
  
    return gtsm

if __name__ == "__main__":
    cfg = load_config('../config.yml') #load config
    
    input_locations = open_input_locations(cfg['input']['paths']['input_locations'],0)
    
    #output_dir = '/Volumes/Naamloos/PhD_Data/GESLA3/preprocessed_dmax/26_03_2024_detrended_deseason_coastal_gesla3_minyr30'
    output_dir = '/Volumes/Naamloos/PhD_Data/GESLA3/preprocessed_dmax/26_03_2024_detrended_deseason_coastal_gtsm_minyr30'
    
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
        
        
    gtsm = save_preprocessed_gtsm(cfg['input']['paths']['gtsm_dmax'],input_locations,
                                         cfg['preprocessing'],cfg['input']['match_dist_limit'],
                                         output_dir=output_dir)
    
    gtsm.to_netcdf(os.path.join(output_dir,'gtsm_preprocessed_dmax_at_gesla3_tgs.nc'))
