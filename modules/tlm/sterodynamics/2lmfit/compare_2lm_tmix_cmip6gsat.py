#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate 2-layer model output (tmix, tdeep, ohc) for each CMIP6 model using
CMIP6 calibration parameters and SSPs ERFs from RCMIP csv

Compare 2-layer model tmix against GSAT from CMIP6 (GSAT obtained from Matt Palmer and Chris Jones, only available for ssp126, 245 and 585)
Store RMSE over these scenarios combined in a NetCDF. 

Input parameters: 
scenarios           = target SSPs
path                = path with 2LM code, input files, etc.
gsat_file           = path to file with CMIP6 GSAT
targ_years          = target years for comparison 2lm tmix and cmip6 gsat
ref_years           = reference years for output

Output:
rmses = root mean square errors across the three scenarios for each individual CMIP6 model with GSAT available (K) [nmodels]

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
gsat_file = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/tas_CMIP6_n17_1986_2005ref_1850_2100_am.nc'
targ_years = np.arange(2015,2101) #target years
ref_years = np.arange(1995,2015) #reference years

erfs,erfyears = fetch_erfs_from_rcmip(csv_path,scenarios)

plt.figure()
for s,scen in enumerate(scenarios):
    plt.plot(erfyears,erfs[s,:],label=scen)
plt.ylabel('ERF [w/m2]')
plt.xlabel('Years')
plt.legend()

#retrieve CMIP6 calibration parameters
scmpy2L_calib = open(path+'scmpy2L_calib_n=44_eps=fit_v20200702.txt')
x = scmpy2L_calib.readlines()[1:]
cmip6_f4x = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[1]) for i in x])))
cmip6_lambda = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[2]) for i in x])))
cmip6_cmix = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[3]) for i in x])))
cmip6_cdeep = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[4]) for i in x])))
cmip6_gamma = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[5]) for i in x])))
cmip6_epsilon = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[6]) for i in x])))

#load CMIP6 gsat
cmip6_gsat = xr.open_dataset(gsat_file,decode_times=True)

#find overlap in models with GSAT and calibration parameters
models = sorted( list(set(list(cmip6_f4x.keys())) & set(cmip6_gsat.model.values)) )

#plotting
fig = plt.figure(figsize=(8,30))
numcols = 4 #number of columns plot

rmses = [] #initialize 
for i,model in enumerate(models):
    
    #add subplot panel
    if i==0: #add subplot to figure handle            
        ax = fig.add_subplot(100+numcols*10+1)
    else: #if not the first model, adjust all previous subplots
        n = len(fig.axes)
        for j in range(n):
            fig.axes[j].change_geometry(int(np.ceil((n+1)/numcols)), numcols, j+1)
            # add the new subplot
            ax = fig.add_subplot(int(np.ceil((n+1)/numcols)), numcols, n+1)
    
    sqerr = 0  
       
    for s,scen in enumerate(['ssp126','ssp245','ssp585']): #for now using just these 3, that's what I have for GSAT
        #define 2lm
        scm2co2 = scmpy2l.ScmDef(
            extforce=erfs[s*2,:], #assume ERF is equal for each model
            exttime=erfyears,
            t2x=None,
            f2x=cmip6_f4x[model]/2,#None ##this should be q2x, but then None throws an error. If you input a value different from the default, the results don't change->only lamg is used as input?
            lamg=cmip6_lambda[model],
            cmix=cmip6_cmix[model],
            cdeep=cmip6_cdeep[model],
            gamma_2l=cmip6_gamma[model],
            eff=cmip6_epsilon[model], tbeg=1850, tend=2100,dt=1)
        
        out2co2 = scm2co2.run()
        
        #reference temperature to reference years and select target years
        gsat_2lm = out2co2.tlev[:,0]-np.mean(out2co2.tlev[:,0][np.where((out2co2.time<=ref_years[-1]) & (out2co2.time>=ref_years[0]))])
        gsat_2lm = gsat_2lm[np.where((out2co2.time<=targ_years[-1]) & (out2co2.time>=targ_years[0]))]
        
        gsat_cmip = cmip6_gsat.sel(model=model).sel(scen=scen).sel(year=targ_years).tas.values-cmip6_gsat.sel(model=model).sel(scen=scen).sel(year=ref_years).tas.mean(dim='year').values
        
        #plot cmip6 and 2lm gsat
        ax.plot(targ_years, gsat_2lm,label='2lm '+scen)     
        ax.plot(targ_years, gsat_cmip,label='cmip6 '+scen) 
        sqerr = sqerr + np.nansum((gsat_cmip-gsat_2lm)**2)
        
    ax.set_title(model+', RMSE='+str('{0:.2f}'.format(np.sqrt(np.mean(sqerr))))+'K')
    ax.set_ylabel('GSAT [K]')    
    ax.set_xlabel('')
   
    rmses = np.append(rmses,np.sqrt(np.mean(sqerr)))     
        
    if i<len(models)-numcols: 
        ax.set_xticklabels('')
    
ax.legend()
plt.show()

#save rmse to netcdf
rmse_da = xr.DataArray(rmses, coords=[models],dims=['model'])
rmse_da.attrs['unit']='K'
rmse_da.attrs['base_period']=str(ref_years[0])+'-'+str(ref_years[-1])
rmse_da.attrs['target_period']=str(targ_years[0])+'-'+str(targ_years[-1])
rmse_da.attrs['long_name']='RMSE_2LMtmix_vs_CMIP6gsat'

ds = rmse_da.to_dataset(name='gsat_rmse')
ds.attrs['info'] = 'RMSE of scmpy 2-layer model mixed layer temperature against GSAT for CMIP6 models obtained from Matt Palmer & Chris Jones across ssp126, ssp245 and ssp585 combined, using CMIP6 calibration parameters with ERF from RCMIP'
ds.to_netcdf(path=path+'scmpy2LM_RCMIP_CMIP6calpm_n'+str(len(models))+'_gsat_rmse.nc',mode='w') #save to file    

#plot CDF
plt.figure(3)
plt.step(np.sort(rmses), np.arange(0,rmses.size)/rmses.size)
plt.axhline(.85,color='black')
plt.ylabel('P [-]') 
plt.xlabel('RMSE [K]')
plt.title('CDF CMIP6-2LM GSAT RMSEs, period: '+str(targ_years[0])+'-'+str(targ_years[-1])+', ref: '+str(ref_years[0])+'-'+str(ref_years[-1]))
plt.text(13,.86,'85%',color='black')
plt.show() 

