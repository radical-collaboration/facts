import xarray as xr
import pandas as pd
import os
import numpy as np
from I_O import open_input_locations, save_ds_to_netcdf, get_refFreqs, lazy_output_to_ds
import argparse
import dask
from esl_analysis import multivariate_normal_gpd_samples_from_covmat, get_return_curve_gpd
from projecting import compute_AFs, compute_AF_timing
from utils import if_scalar_to_list
from dask.distributed import LocalCluster

# Do not warn about chained assignments
pd.options.mode.chained_assignment = None  # default='warn'

""" extremesl_fit.py
to-do: document IO
"""

def compute_projectESL_output(loc,scale,shape,rate,cov,mhhw,f,below_threshold,n_samples,refFreq,input_locations,out_qnts,target_years,target_AFs,target_freqs,z_hist=None):
    if z_hist is not None: #if return curve already provided
        z = z_hist
        
    else: #derive return curves from GPD samples
        if cov is not None: #generate scale and shape samples
            scale_samples,shape_samples = multivariate_normal_gpd_samples_from_covmat(scale,shape,cov,n_samples,0) 
        else: #use best estimates (means) or provided samples
            scale_samples,shape_samples = scale,shape #use central estimate scale & shape
 
        z= get_return_curve_gpd(f,scale_samples,shape_samples,loc,rate,below_threshold,mhhw) #compute return curve samples
       
    #compute quantiles of output:
    if z.ndim > 1:
        z_hist_quantiles = np.quantile(z,out_qnts,axis=-1) #dont have to use nanquantile, because either defined for all samples or for none
    else: 
        z_hist_quantiles = np.repeat(z[None,:],len(out_qnts),axis=0) #central estimate parameters used to compute z, so there is no uncertainty in z_hist --> repeat
        
    if len(target_years)>0: #compute AFs for target years
        af,max_af,z_fut = compute_AFs(f,z,input_locations.sea_level_change.sel(years=target_years),refFreq)
        z_fut_quantiles = np.quantile(z_fut,out_qnts,axis=1)    
        af_quantiles = np.quantile(af,out_qnts,axis=0)
    else:
        af_quantiles = np.nan
        max_af = np.nan
        z_fut_quantiles = np.nan
        
    if len(target_AFs)>0: #compute timing of target AFs
        af_timing = []
        for target_af in np.array(target_AFs): #loop over target AFs to compute timing
            af_timing.append(compute_AF_timing(f,z,input_locations.sea_level_change,refFreq,target_af))
            
        af_timing_quantiles = np.quantile(np.vstack((af_timing)),out_qnts,axis=-1)
        af_timing_quantiles = np.round(np.where(af_timing_quantiles>=input_locations.years[-1].values,9999,af_timing_quantiles))
    else:
        af_timing_quantiles = np.nan
     
    if len(target_freqs)>0: #compute timing of target AFs
        f_timing = []
        for target_freq in np.array(target_freqs): #loop over target AFs to compute timing
            f_timing.append(compute_AF_timing(f,z,input_locations.sea_level_change,refFreq,np.round(target_freq/refFreq).astype('int')))
            
        f_timing_quantiles = np.quantile(np.vstack((f_timing)),out_qnts,axis=-1)
        f_timing_quantiles = np.round(np.where(f_timing_quantiles>=input_locations.years[-1].values,9999,f_timing_quantiles))
    else:
        f_timing_quantiles = np.nan
        
    return z_hist_quantiles,z_fut_quantiles,af_quantiles,max_af,af_timing_quantiles,f_timing_quantiles


def project_ESLs_lazily(esl_statistics,f,below_threshold,n_samples,refFreqs,input_locations,out_qnts,target_years,target_AFs,target_freqs):
    target_years = if_scalar_to_list(target_years)
    target_AFs = if_scalar_to_list(target_AFs)
    target_freqs = if_scalar_to_list(target_freqs)
    
    locs = scales = shapes = rates = covs = mhhws = z_hists = [None for k in esl_statistics.locations] #initialize input parameters for each location with None
    
    if 'z_hist' in esl_statistics: #if return curves provided, use these as input
        z_hists = esl_statistics['z_hist'].values
    
    else: #else use GPD parameters as input
        locs    = esl_statistics['loc'].values
        scales  = esl_statistics['scale'].values
        shapes  = esl_statistics['shape'].values
        rates   = esl_statistics['avg_extr_pyear'].values
        
        if 'cov' in esl_statistics: # if covariance matrix is provided, use it
            covs = [k.tolist() for k in esl_statistics['cov'].values]
            
        elif 'scale_samples' in esl_statistics and 'shape_samples' in esl_statistics: #if scale/shape samples provided directly
            scales = esl_statistics['scale_samples'].values
            shapes = esl_statistics['shape_samples'].values
         
        if 'mhhw' in esl_statistics:
            mhhws = esl_statistics['mhhw'].values
        else:
            if below_threshold == 'mhhw':
                print('Warning: cannot compute return heights below threshold using "mhhw" because MHHW value is not available from ESL data. Continuing without modeling below threshold.')
                below_threshold = None
    lazy_results=[]
    
    for l,location in enumerate(esl_statistics.locations):
        lazy_result = dask.delayed(
                            compute_projectESL_output)(locs[l],scales[l],shapes[l],rates[l],covs[l],mhhws[l],f,below_threshold,n_samples,
                                                       refFreqs[l],input_locations.sel(locations=location),out_qnts,
                                                       target_years,target_AFs,target_freqs,z_hists[l])
        lazy_results = np.append(lazy_results,lazy_result)
    
    return lazy_results

if __name__ == "__main__":
    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the projection stage for the extreme sea-level workflow",
                    epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
    
    # Define the command line arguments to be expected
    parser.add_argument('--nsamps', help="Number of samples to draw [default = 20000]", type=int, default=20000)
    parser.add_argument('--out_qnts', help="Which quantiles to store output for", type=float, default=[.05,.5,.95]) #the following for inputs are arrays, probably not possible to parse?
    parser.add_argument('--target_years', help="")
    parser.add_argument('--target_AFs', help="" )
    parser.add_argument('--target_freqs', help="")
    parser.add_argument('--total_localsl_file', help="Total localized sea-level projection file. Site lats/lons are taken from this file and mapped to the GESLA database", default="total-workflow_localsl.nc")
    parser.add_argument('--refFreq_data', help="Which protection level frequencies data to use", default="flopros")
    parser.add_argument('--refFreq_path', help="Directory containing requested protection level frequencies data", default=os.path.join(os.path.dirname(__file__), ""))
    parser.add_argument('--use_central_esl_estimates_only', help="Whether to use uncertainty of ESL estimates for projections", default=False)
    parser.add_argument('--below_threshold', help="How to model events below the GPD threshold", default='mhhw')
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
    
    # Parse the arguments
    args = parser.parse_args()
    
    sl_fn               = args.total_localsl_file
    n_samples           = args.nsamps
    out_qnts            = args.out_qnts
    path_to_refFreqs    = args.refFreq_path
    refFreq_data        = args.refFreq_data
    use_central_esl_estimates_only        = args.use_central_esl_estimates_only
    below_threshold = args.below_threshold
    target_AFs = args.target_AFs
    target_years = args.target_years
    target_freqs = args.yarget_freqs
    
    cluster = LocalCluster()          # Start a fully-featured local Dask cluster
    client = cluster.get_client()
    
    f= 10**np.linspace(-6,2,num=1001) #input frequencies to compute return heights for
    f=np.append(f,np.arange(101,183))    
    
    input_locations = open_input_locations(sl_fn,n_samples)
    esl_statistics = xr.open_dataset('esl_statistics.nc')
    
    ### projecting stage
    refFreqs = get_refFreqs(refFreq_data,input_locations,esl_statistics,path_to_refFreqs) #grab reference frequencies for AFs at each site

    if use_central_esl_estimates_only:
        esl_statistics=esl_statistics.drop(['cov','scale_samples','shape_samples'],errors='ignore')
        
    output=dask.compute(*project_ESLs_lazily(esl_statistics,f,below_threshold,n_samples,refFreqs,
                                         input_locations,out_qnts,target_years,target_AFs,target_freqs)) #compute output for each location in parallel

    output_ds = lazy_output_to_ds(output,f,out_qnts,esl_statistics,target_years,target_AFs,target_freqs) #convert output to xarray dataset
    output_ds['refFreq'] = (['locations'],refFreqs)

    output_ds.to_netcdf('projectESL_output.nc','w') #store
    exit()