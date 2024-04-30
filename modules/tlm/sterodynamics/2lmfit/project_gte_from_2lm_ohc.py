#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate GTE samples from OHC samples using expansion coefficients fitted
to 2-layer model output and CMIP6 GTE. The expansion coefficient distribution 
is clipped based on the GSAT rmse, and then assumed normally distributed.

Created on Mon Sep 14 16:08:01 2020

@author: thermans
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm

scenarios =['ssp119','ssp126','ssp245','ssp370','ssp585'] #target scenarios
path = '/Users/thermans/Documents/IPCC_AR6/2LM/' #path

pb_clip = .85 #some metric to cut off the GSAT rmse CDF at

seed=1234
np.random.seed(seed)

#load OHC samples
ohc_ds = xr.open_dataset(path+'scmpy2LM_RCMIP_CMIP6calpm_ens_t_ohc.nc')
ohc_samps = ohc_ds.ohc.values

#load expansion coefficients
expcoefs = xr.open_dataset(path+'scmpy2LM_RCMIP_CMIP6calpm_n18_expcoefs.nc').expcoefs

#load GSAT RMSEs
rmses = xr.open_dataset(path+'scmpy2LM_RCMIP_CMIP6calpm_n17_gsat_rmse.nc').gsat_rmse

#clip gsat rmse cdf
rmses_sorted = np.sort(rmses)
sort_idx = np.argsort(rmses)
pbs = np.arange(0,rmses_sorted.size)/rmses_sorted.size

expcoefs_clipped=expcoefs.drop(rmses.model[sort_idx.values][pbs>pb_clip].values,dim='model') #caveat: gsat rmses and expcoefs do not have exactly the same models, because of data availability

#extract normal distribution modes
mean_expcoefs = expcoefs_clipped.mean(dim='model')
std_expcoefs = expcoefs_clipped.std(dim='model')

#generate samples assuming normal distribution
expcoef_samps = np.random.normal(loc=mean_expcoefs, scale=std_expcoefs, size=ohc_samps.shape)
gte_samps=ohc_samps*expcoef_samps #the GTE output


#plot some diagnostics
#gte
plt.figure()
colors = ['C1','C0','C2','C3','C4']
plt.grid()

for s,scen in enumerate(scenarios):
    #plot median and 5-95%    
    plt.plot(ohc_ds.years,np.median(gte_samps[s,:,:],axis=0),label=scenarios[s],color=colors[s])
    plt.fill_between(ohc_ds.years,np.quantile(gte_samps[s,:,:],.05,axis=0),np.quantile(gte_samps[s,:,:],.95,axis=0),color=colors[s],alpha=0.3)

plt.legend()
plt.xlabel('years')
plt.ylabel('GTE [m]')

#histograms of expansion coefficients
(mu1, sigma1) = norm.fit(expcoefs)
(mu2, sigma2) = norm.fit(expcoefs_clipped)

fig,axs = plt.subplots(1,2)
axs[0].hist(expcoefs.values,bins=8)
axs[0].set_ylabel('n'); axs[0].set_xlabel('Expansion Coefficient [m/YJ]')
axs[0].set_title('all models')
axs[0].set_xlim(.065,.14)
axs[0].set_ylim(0,7)
axs[0].text(.08,6.5,'mean='+str('{0:.3f}'.format(mu1))+', std='+str('{0:.3f}'.format(sigma1)))


axs[1].hist(expcoefs_clipped.values,bins=6)
axs[1].set_ylabel('n'); axs[1].set_xlabel('Expansion Coefficient [m/YJ]')
axs[1].set_title('clipped distribution')
axs[1].set_xlim(.065,.14)
axs[1].set_ylim(0,7)
axs[1].text(.08,6.5,'mean='+str('{0:.3f}'.format(mu2))+', std='+str('{0:.3f}'.format(sigma2)))
