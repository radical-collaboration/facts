import xarray as xr
import dask
import pandas as pd
import os
import numpy as np
import argparse
from I_O import load_config, open_input_locations, esl_statistics_dict_to_ds, save_ds_to_netcdf, open_gpd_parameters,get_coast_rp_return_curves
from esl_analysis import ESL_stats_from_raw_GESLA, ESL_stats_from_gtsm_dmax, multivariate_normal_gpd_samples_from_covmat, get_return_curve_gpd
import tarfile
import fnmatch

# Do not warn about chained assignments
pd.options.mode.chained_assignment = None  # default='warn'

""" extremesl_fit.py
to-do: document IO
"""

def get_ESL_statistics(esl_data,path_to_data,input_locations,match_dist_limit,preproc_settings=None,n_samples=None,f=None):
    print('Extracting ESL information for queried input locations...')
    esl_statistics = {} #initialize dictionary to hold ESL information
    esl_data = esl_data.lower()
    
    if esl_data in ['gesla2','gesla3']: #if using raw data from GESLA
        esl_statistics = ESL_stats_from_raw_GESLA(esl_data,path_to_data,input_locations,preproc_settings,match_dist_limit) #dictionary output
        esl_statistics = esl_statistics_dict_to_ds(input_locations,esl_statistics)
        
    elif esl_data.lower() in ['gtsm_dmax']:
        esl_statistics = ESL_stats_from_gtsm_dmax(path_to_data,input_locations,preproc_settings,match_dist_limit)
            
    elif esl_data.lower() in ['hermans2023','kirezci2020','vousdoukas2018','gtsm_dmax_gpd']:
        esl_statistics = open_gpd_parameters(esl_data,path_to_data,input_locations,n_samples,match_dist_limit)
            
    elif esl_data.lower() == 'coast-rp': # Option C: read in pre-defined return curves         
        esl_statistics = get_coast_rp_return_curves(path_to_data,input_locations,f,match_dist_limit)
    else:
        raise Exception('ESL input data type not recognized.')
    return esl_statistics

if __name__ == "__main__":
    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the fitting stage for the extreme sea-level workflow",
                    epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
    
    # Define the command line arguments to be expected
    parser.add_argument('--minYears', help="Minimum number of years available [default=20]", type=int, default=20)
    parser.add_argument('--resample_freq', help="Frequency to resample the raw data to prior to ESL analysis.", default="D_max")
    parser.add_argument('--deseasonalize', help="Boolean flag to indicate whether to remove mean seasonal cycle prior to ESL analysis.",type=int, default=1)
    parser.add_argument('--detrend', help="Boolean flag to indicate whether to remove linear trend prior to ESL analysis.",type=int, default=1)
    parser.add_argument('--subtract_amean', help="Boolean flag to indicate whether to remove annual means prior to ESL analysis.",type=int, default=1)
    parser.add_argument('--match_lim', help="Radius around requested locations to find a matching tide gauge in GESLA database", type=float, default=10)
    parser.add_argument('--gpd_pot_threshold', help="Percentile for GPD analysis [default=99]", type=float, default=99)
    parser.add_argument('--decluster_window', help="Maximum number of days that define a cluster for extreme events [default=3]", type=int, default=3)
    parser.add_argument('--decluster_method', help="Method to use for declustering peaks.", default="rolling_max")
    parser.add_argument('--nsamps', help="Number of samples to draw [default = 20000]", type=int, default=2000)
    parser.add_argument('--total_localsl_file', help="Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database", default="total-workflow_localsl.nc")
    parser.add_argument('--esl_data', help="Type of data used for the ESL analysis.", default="gesla3")
    parser.add_argument('--esl_data_path', help="Directory containing requested ESL data", default=os.path.join(os.path.dirname(__file__), "gesla3_data"))
	
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
    
    # Parse the arguments
    args = parser.parse_args()
    
    #generate dictionary with preprocessing settings
    preproc_settings = {}
    preproc_settings['min_yrs']             = args.minYears
    preproc_settings['store_esls']          = False
    preproc_settings['resample_freq']       = args.resample_freq
    preproc_settings['deseasonalize']       = args.deseasonalize
    preproc_settings['detrend']             = args.detrend
    preproc_settings['subtract_amean']      = args.subtract_amean
    preproc_settings['ref_to_msl']          = 0
    preproc_settings['declus_method']       = args.decluster_method
    preproc_settings['declus_window']       = args.decluster_window
    preproc_settings['extremes_threshold']  = args.gpd_pot_threshold                                             
    
    esl_data        = args.esl_data
    esl_data_path   = args.esl_data_path
    sl_fn           = args.total_localsl_file
    n_samples       = args.nsamps
    
    f= 10**np.linspace(-6,2,num=1001) #input frequencies to compute return heights for
    f=np.append(f,np.arange(101,183))    
    
    input_locations = open_input_locations(sl_fn,n_samples)
    extremesl_fit = get_ESL_statistics(esl_data,esl_data_path,input_locations,args.match_lim,preproc_settings,n_samples,f)
    extremesl_fit.to_netcdf(os.path.join(os.path.dirname(__file__),'{}_esl_statistics.nc'.format(args.pipeline_id)),mode='w')
    exit()