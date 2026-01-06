#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" regress_cmip6GTE_to_2lmOHC.py
This script runs the regression task for the 2-layer model global thermal 
expansion analysis. This task finds expansion coefficients by linear regression
of dedrifted CMIP6 'zostoga' output against the ocean heat content output of 
Chris Smith's simple two-layer model code. The input to the 2-layer model are
the calibration parameters (provided by CH7) of the CMIP6 models, and
model-independent effective radiative forcing from RCMIP (https://www.rcmip.org/).

The regression is 1D in the sence that we use the total OHC (mix+deep), and
scenario-independent in the sense that the data across scenarios is pooled before 
regression. 

Input parameters: 
scenarios           = target SSPs
path                = path with 2LM code, input files, etc.
gte_file            = path to file with CMIP6 dedrifted GTE
targ_years          = target years for regression
ref_years           = reference years for regression

Output:
cmip6_expcoefs      = CMIP6 model-dependent expansion coefficients (m/YJ) [nmodels]

Created on Fri Sept 25 09:44:27 2020
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
from twolayermodel.scmpy import scmpy2l   
import xarray as xr
from scipy.optimize import curve_fit

plt.close('all')

def linfunc(x,a): #fitting function without intercept
    return a*x

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
scenarios =['ssp119','ssp126','ssp245','ssp370','ssp585'] #target scenarios
path = '/Users/thermans/Documents/IPCC_AR6/2LM/' #path
gte_file = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/zostoga_CMIP6_AR6_1986_2005ref_ldedr_1850_2100_am.nc' 
targ_years = np.arange(2015,2101) #target years for regression
ref_years = np.arange(1995,2015) #reference years for regression

erfs,erfyears = fetch_erfs_from_rcmip(path,scenarios) #get ERF timeseries

#retrieve CMIP6 calibration parameters
scmpy2L_calib = open(path+'scmpy2L_calib_n=44_eps=fit_v20200702.txt')
x = scmpy2L_calib.readlines()[1:]
cmip6_f4x = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[1]) for i in x])))
cmip6_lambda = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[2]) for i in x])))
cmip6_cmix = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[3]) for i in x])))
cmip6_cdeep = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[4]) for i in x])))
cmip6_gamma = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[5]) for i in x])))
cmip6_epsilon = dict(zip([i.split()[0] for i in x], np.array([float(i.split()[6]) for i in x])))


#get dedrifted CMIP6 'zostoga'
cmip6_gte = xr.open_dataset(gte_file,decode_times=True) 

#find overlap in models with GTE and calibration parameters
models = sorted( list(set(list(cmip6_f4x.keys())) & set(cmip6_gte.model.values)) )

cmip6_expcoefs = np.empty((len(models)))  #initialize 
for m,model in enumerate(models):
    #initialize empty arrays to store ohc and gte across all scenarios
    model_ohc=np.empty((0))      
    model_gte=np.empty((0))  
             
    for s,scen in enumerate(['ssp119','ssp126','ssp245','ssp370','ssp585']): #loop over scenarios
        scm2co2 = scmpy2l.ScmDef( #define 2lm
            extforce=erfs[s,:],
            exttime=erfyears,
            t2x=None,
            f2x=None,
            lamg=cmip6_lambda[model], #feedback parameter
            cmix=cmip6_cmix[model], #heat capacity surface
            cdeep=cmip6_cdeep[model], #heat capacity deep
            gamma_2l=cmip6_gamma[model], #heat exchange coefficient
            eff=cmip6_epsilon[model], tbeg=1850, tend=2100,dt=1) #efficacy
        
        out2co2 = scm2co2.run() #run 2lm
        
        #convert OHC units from 10**22 to 10**24 (YJ), and reference temperature and OHC output to reference years
        out2co2.ohc = 0.01*(out2co2.ohc - np.mean(out2co2.ohc[np.where((out2co2.time<=ref_years[-1]) & (out2co2.time>=ref_years[0]))]))
        out2co2.ohc = out2co2.ohc[np.where((out2co2.time<=targ_years[-1]) & (out2co2.time>=targ_years[0]))] #grab target years
        
        #store OHC and GTE in model array, if scenario exists for model
        if np.isnan(cmip6_gte.sel(model=model).sel(scen=scen).zostoga.values).any():
            continue
        else:
            #group results for all scenarios, relative to reference years
            model_ohc = np.append(model_ohc,out2co2.ohc) 
            model_gte = np.append(model_gte,cmip6_gte.sel(model=model).sel(scen=scen).sel(year=targ_years).zostoga.values
                                  -cmip6_gte.sel(model=model).sel(scen=scen).sel(year=ref_years).zostoga.mean(dim='year').values) #select target years and subtract base

    #linear regression without intercept for all scenarios combined        
    slope, cov = curve_fit(linfunc, model_ohc, model_gte)
    cmip6_expcoefs[m] = slope #m/YJ

#save expansion coefficients to netcdf
expcoefs_da = xr.DataArray(cmip6_expcoefs, coords=[models],dims=['model'])
expcoefs_da.attrs['unit']='m/YJ'
expcoefs_da.attrs['base_period']=str(ref_years[0])+'-'+str(ref_years[-1])
expcoefs_da.attrs['target_period']=str(targ_years[0])+'-'+str(targ_years[-1])
expcoefs_da.attrs['long_name']='expansion_coefficients'

ds = expcoefs_da.to_dataset(name='expcoefs')
ds.attrs['info'] = 'Expansion coefficients derived by regressing scmpy 2-layer model OHC against linearly dedrifted GTE for CMIP6 models, using CMIP6 calibration parameters with ERF from RCMIP'
ds.to_netcdf(path=path+'scmpy2LM_RCMIP_CMIP6calpm_n'+str(len(models))+'_expcoefs.nc',mode='w') #save to file