#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:40:17 2023

@author: timhermans
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

a = xr.open_dataset('/Volumes/Naamloos/PhD_Data/GESLA3/timing_projections/esl_AF_timing_ar6_2150wfs_yrRef2022_fRefFLOPROS_muATS_extrapSweet.nc')
b = xr.open_dataset('/Volumes/Naamloos/PhD_Data/COAST_RP/COAST-RP.nc')

plt.figure()
plt.plot(a.freq,a.rc_ce.isel(site=1))

plt.xscale('log')
plt.ylim([-1,6])
plt.xlim([1e-6,365.25/2])
plt.yticks(np.arange(-1,6.5,.5))
plt.grid()
