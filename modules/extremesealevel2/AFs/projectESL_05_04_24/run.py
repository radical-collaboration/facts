#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:51:03 2023

Execution of projectESL according to user configuration in 'config.yml'.

@author: timhermans
"""
import xarray as xr
import dask
import os
import numpy as np
from I_O import load_config, get_refFreqs, open_input_locations, lazy_output_to_ds, esl_statistics_dict_to_ds, save_ds_to_netcdf
from I_O import open_gpd_parameters,get_coast_rp_return_curves
from esl_analysis import ESL_stats_from_raw_GESLA, ESL_stats_from_gtsm_dmax, multivariate_normal_gpd_samples_from_covmat, get_return_curve_gpd
from projecting import compute_AFs, compute_AF_timing
from utils import if_scalar_to_list

def get_ESL_statistics(esl_data,path_to_data,input_locations,match_dist_limit,preproc_settings=None,n_samples=None,f=None,output_dir=None):
    print('Extracting ESL information for queried input locations...')
    esl_statistics = {} #initialize dictionary to hold ESL information
    esl_data = esl_data.lower()
    
    if esl_data in ['gesla2','gesla3']: #if using raw data from GESLA
        esl_statistics = ESL_stats_from_raw_GESLA(esl_data,path_to_data,input_locations,preproc_settings,match_dist_limit,output_dir) #dictionary output
        esl_statistics = esl_statistics_dict_to_ds(input_locations,esl_statistics)
        
    elif cfg['input']['esl_data'].lower() in ['gtsm_dmax']:
        esl_statistics = ESL_stats_from_gtsm_dmax(path_to_data,input_locations,preproc_settings,match_dist_limit,output_dir)
            
    elif cfg['input']['esl_data'].lower() in ['hermans2023','kirezci2020','vousdoukas2018','gtsm_dmax_gpd']:
        esl_statistics = open_gpd_parameters(esl_data,path_to_data,input_locations,n_samples,match_dist_limit)
            
    elif cfg['input']['esl_data'].lower() == 'coast-rp': # Option C: read in pre-defined return curves         
        esl_statistics = get_coast_rp_return_curves(path_to_data,input_locations,f,match_dist_limit)
    else:
        raise Exception('ESL input data type not recognized.')
    return esl_statistics


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
#### 

if __name__ == "__main__":
    from dask.distributed import LocalCluster
    cluster = LocalCluster()          # Start a fully-featured local Dask cluster
    client = cluster.get_client()
    
    cfg = load_config('../config.yml') #load config
    
    output_dir = os.path.join(cfg['general']['output_dir'],cfg['general']['run_name'])
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
        
    #assign a few settings to local variables
    n_samples               = cfg['general']['n_samples'] #number of samples for uncertainty propagation
    esl_data                = cfg['input']['esl_data'] #type of ESL data to load/analyze
    out_qnts                = np.array(cfg['general']['output_quantiles'].split(',')).astype('float') #output quantiles to evaluate results at
    refFreq_data            = cfg['projecting']['refFreqs'] #type of refFreq data to load (or scalar to use at all locations)
    below_threshold         = cfg['projecting']['below_threshold']
    match_dist_limit        = cfg['input']['match_dist_limit']
    if below_threshold == 'None':
        below_threshold = None
        
    if refFreq_data in ['diva','flopros']:
        path_to_refFreqs    = cfg['input']['paths'][refFreq_data] #path to refFreqs if applicable
    
    if cfg['projecting']['target_years'] is not None:
        target_years = np.array(str(cfg['projecting']['target_years']).split(',')).astype('int')
    else:
        target_years = []
      
    if cfg['projecting']['target_AFs'] is not None:
        target_AFs = np.array(str(cfg['projecting']['target_AFs']).split(',')).astype('int')
    else:
        target_AFs = []
        
    if cfg['projecting']['target_freqs'] is not None:
        target_freqs = np.array(str(cfg['projecting']['target_freqs']).split(',')).astype('float')  
    else:
        target_freqs = []
    
    f= 10**np.linspace(-6,2,num=1001) #input frequencies to compute return heights for
    f=np.append(f,np.arange(101,183))     
        
    #open input locations
    input_locations = open_input_locations(cfg['input']['paths']['input_locations'],n_samples)
    #input_locations = input_locations.isel(locations = np.where(np.isfinite(input_locations.sea_level_change.isel(years=0,samples=0)))[0])#temporary, fix nans in nearest interpolation from full gridded samples
    
    ### fitting stage
    esl_statistics = get_ESL_statistics(esl_data,cfg['input']['paths'][esl_data],input_locations,match_dist_limit,cfg['preprocessing'],n_samples,f,output_dir)  #get ESL information at input locations
    esl_statistics.attrs['cfg'] = str(cfg)
    save_ds_to_netcdf(output_dir,esl_statistics,'esl_statistics.nc') #store
    
    ### projecting stage
    refFreqs = get_refFreqs(refFreq_data,input_locations,esl_statistics,path_to_refFreqs) #grab reference frequencies for AFs at each site
    
    if cfg['general']['use_central_esl_estimates_only']:
        esl_statistics=esl_statistics.drop(['cov','scale_samples','shape_samples'],errors='ignore')
        
    output=dask.compute(*project_ESLs_lazily(esl_statistics,f,below_threshold,n_samples,refFreqs,
                                         input_locations,out_qnts,target_years,target_AFs,target_freqs)) #compute output for each location in parallel
    
    output_ds = lazy_output_to_ds(output,f,out_qnts,esl_statistics,target_years,target_AFs,target_freqs) #convert output to xarray dataset
    output_ds['refFreq'] = (['locations'],refFreqs)
    output_ds.attrs['cfg'] = str(cfg)
    
    save_ds_to_netcdf(output_dir,output_ds,'projectESL_output.nc') #store

    cluster.close()
    '''
    import matplotlib.pyplot as plt
  
    plt.figure()
    ax = plt.subplot(211)
    output_ds.isel(locations=0).z_hist.plot.line(x='f',ax=ax)
    output_ds.isel(locations=0).z_fut.sel(target_year=2100).plot.line(x='f',ax=ax,color='blue')
    ax.set_xscale('log')
    ax.set_xlabel('Return frequency [1/yr]')
    ax.set_ylabel('Return height [m]')
    ax.axvline(x=output_ds.isel(locations=0).refFreq,color='grey',linestyle='dashed')
    ax.axvline(x=output_ds.isel(locations=0).refFreq * output_ds.isel(locations=0,target_year=-1).AF.values[0],color='grey',linestyle='dotted')
    ax.axvline(x=output_ds.isel(locations=0).refFreq * output_ds.isel(locations=0,target_year=-1).AF.values[1],color='grey',linestyle='dotted')
    ax.axvline(x=output_ds.isel(locations=0).refFreq * output_ds.isel(locations=0,target_year=-1).AF.values[2],color='grey',linestyle='dotted')
    
    ax = plt.subplot(212)
    input_locations.sea_level_change.load().isel(locations=0).quantile(out_qnts,dim='samples').plot.line(x='years',ax=ax)
    
    plt.figure()
    plt.plot(output_ds.isel(locations=0).z_hist.isel(qnt=1),f,color='black')
    plt.plot(output_ds.isel(locations=0).z_hist.isel(qnt=0),f,color='black',linestyle='dashed')
    plt.plot(output_ds.isel(locations=0).z_hist.isel(qnt=-1),f,color='black',linestyle='dashed')
    
    plt.plot(output_ds.isel(locations=0).z_fut.isel(qnt=1).sel(target_year=2100),f,color='red')
    plt.plot(output_ds.isel(locations=0).z_fut.isel(qnt=0).sel(target_year=2100),f,color='red',linestyle='dashed')
    plt.plot(output_ds.isel(locations=0).z_fut.isel(qnt=-1).sel(target_year=2100),f,color='red',linestyle='dashed')
    plt.axhline(y=1e-2)
    plt.axhline(y=1e-3)
    plt.yscale('log')
    plt.ylim([1e-6,1e2])
    plt.xlim([0,6])
    '''