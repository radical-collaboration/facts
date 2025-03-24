#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate 2-layer model ensemble output (tmix, tdeep, ohc) by random sampling from 
CMIP6 calibration parameters, using SSPs ERFs from RCMIP csv

Input parameters: 
scenarios           = target SSPs
path                = path with 2LM code, input files, etc.
targ_years          = target years for generation of 2lm output
ref_years           = reference years for output
nsamps              = number of samples
seed                = random number generator seed

Output:
tmix_ens     = mixed layer temperature (K) [nscen,nsamps,nyears]
tdeep_ens    = deep ocean temperature (K) [nscen,nsamps,nyears]
ohc_ens      = ocean heat content (YJ) [nscen,nsamps,nyears]

Created on Mon Aug 10 09:48:55 2020
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
from twolayermodel.scmpy import scmpy2l   
import xarray as xr

plt.close('all')

def fetch_erfs_from_rcmip(csv_path,scenarios):
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
    ssp_idx = {'ssp119':212,'ssp126':231,'ssp245':308,'ssp370':59,'ssp585':404}
    
    with open(path+'rfmip-radiative-forcing-annual-means-v4-0-0.csv') as csv_file:
        csv_reader = csv.reader(csv_file)
        rows = list(csv_reader)
    
    erfs = np.empty((len(scenarios),len(erfyears)))
    
    for s,scen in enumerate(scenarios):
        try:    
            erfs[s,:] = rows[ssp_idx[scen]][7:]
        except:
            continue
        
    return erfs, erfyears

#user settings
scenarios =['ssp119','ssp126','ssp245','ssp370','ssp585']
path = '/Users/thermans/Documents/IPCC_AR6/2LM/'
csv_path = path + 'rfmip-radiative-forcing-annual-means-v4-0-0.csv'

targ_years = np.arange(1850,2301) #target years
ref_years = np.arange(1995,2015) #reference years

erfs,erfyears = fetch_erfs_from_rcmip(csv_path,scenarios)

nsamps = 5000
seed=1234
np.random.seed(seed)

#retrieve CMIP6 calibration parameters
scmpy2L_calib = open(path+'scmpy2L_calib_n=44_eps=fit_v20200702.txt')
x = scmpy2L_calib.readlines()[1:]
cmip6_f4x = np.array([float(i.split()[1]) for i in x])
cmip6_lambda = np.array([float(i.split()[2]) for i in x])
cmip6_cmix = np.array([float(i.split()[3]) for i in x])
cmip6_cdeep = np.array([float(i.split()[4]) for i in x])
cmip6_gamma = np.array([float(i.split()[5]) for i in x])
cmip6_epsilon = np.array([float(i.split()[6]) for i in x])

tmix_ens = np.empty((len(scenarios),nsamps,len(targ_years)))  #initialize 
tdeep_ens = np.empty((len(scenarios),nsamps,len(targ_years)))
ohc_ens = np.empty((len(scenarios),nsamps,len(targ_years)))

#plotting
plt.figure()
colors = ['C1','C0','C2','C3','C4']
plt.grid()
for s,scen in enumerate(scenarios):
    for sp,samp in enumerate(np.arange(nsamps)):
        idx = np.random.randint(0,len(cmip6_lambda)) #index to sample randomly from calibration parameter sets
        #define 2lm
        scm2co2 = scmpy2l.ScmDef(
            extforce=erfs[s,:], #grab ERF for scenario, assume ERF is model independent
            exttime=erfyears, #years corresponding to ERF input
            t2x=None,
            f2x=None, #this should be q2x, but then None throws an error. If you input a value different from the default, the results don't change->only lamg is used as input?
            lamg=cmip6_lambda[idx], #calibration parameters
            cmix=cmip6_cmix[idx],
            cdeep=cmip6_cdeep[idx],
            gamma_2l=cmip6_gamma[idx],
            eff=cmip6_epsilon[idx], tbeg=targ_years[0], tend=targ_years[-1],dt=1)
        
        out2co2 = scm2co2.run()
        
        #convert OHC units from 10**22 to 10**24 (YJ), and reference temperature and OHC output to reference years
        tmix = out2co2.tlev[:,0]-np.mean(out2co2.tlev[:,0][np.where((out2co2.time<=ref_years[-1]) & (out2co2.time>=ref_years[0]))]) #[K]
        tdeep = out2co2.tlev[:,1]-np.mean(out2co2.tlev[:,1][np.where((out2co2.time<=ref_years[-1]) & (out2co2.time>=ref_years[0]))]) #[K]
        ohc = 0.01*(out2co2.ohc-np.mean(out2co2.ohc[np.where((out2co2.time<=ref_years[-1]) & (out2co2.time>=ref_years[0]))])) #[YJ, 10**24]
    
        tmix_ens[s,sp,:] = tmix
        tdeep_ens[s,sp,:] = tdeep
        ohc_ens[s,sp,:] = ohc
        
    #plot median and 5-95%    
    plt.plot(targ_years,np.median(tmix_ens[s,:,:],axis=0),label=scenarios[s],color=colors[s])
    plt.fill_between(targ_years,np.quantile(tmix_ens[s,:,:],.05,axis=0),np.quantile(tmix_ens[s,:,:],.95,axis=0),color=colors[s],alpha=0.3)

plt.legend()
plt.xlabel('years')
plt.ylabel('T [K]')
    
#save to netcdf
tmix_da = xr.DataArray(tmix_ens, coords=[scenarios,np.arange(nsamps),targ_years],dims=['scenario','samps','years'])
tmix_da.attrs['unit']='K'
tmix_da.attrs['base_period']=str(ref_years[0])+'-'+str(ref_years[-1])
tmix_da.attrs['long_name']='mixed_layer_temperature'

tdeep_da = xr.DataArray(tdeep_ens,  coords=[scenarios,np.arange(nsamps),targ_years],dims=['scenario','samps','years'])
tdeep_da.attrs['unit']='K'
tdeep_da.attrs['base_period']=str(ref_years[0])+'-'+str(ref_years[-1])
tdeep_da.attrs['long_name']='deep_ocean_temperature'

ohc_da = xr.DataArray(ohc_ens, coords=[scenarios,np.arange(nsamps),targ_years],dims=['scenario','samps','years'])
ohc_da.attrs['unit']='YJ'
ohc_da.attrs['base_period']=str(ref_years[0])+'-'+str(ref_years[-1])
ohc_da.attrs['long_name']='ocean_heat_content'

ds = tmix_da.to_dataset(name='tmix')
ds['tdeep'] = tdeep_da
ds['ohc'] = ohc_da
ds.attrs['info'] = 'scmpy 2-layer model temperature and ocean heat content output based on random sampling from CMIP6 calibration parameters with ERF from RCMIP'
ds.to_netcdf(path=path+'scmpy2LM_RCMIP_CMIP6calpm_ens_t_ohc.nc',mode='w') #save to file    


