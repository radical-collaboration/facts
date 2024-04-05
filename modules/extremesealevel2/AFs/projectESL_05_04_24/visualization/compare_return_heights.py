#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 08:44:41 2024

compare return heights at the catchments of Jack Heslop

@author: timhermans
"""

import geopandas as geopd
import xarray as xr
import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

plt.close('all')

crp = xr.open_dataset('//Users/timhermans/Documents/GitHub/projectESL/output/COAST_RP/COAST_RP/projectESL_output.nc')

gtsm_ref = xr.open_dataset('/Users/timhermans/Documents/GitHub/projectESL/output/COAST_RP/Muis2023_dmax_from_hourly_refSettings/projectESL_output.nc')
gtsm_99 = xr.open_dataset('/Users/timhermans/Documents/GitHub/projectESL/output/COAST_RP/Muis2023_dmax_from_hourly_99/projectESL_output.nc')

gesla3_ref = xr.open_dataset('/Users/timhermans/Documents/GitHub/projectESL/output/COAST_RP/GESLA3_refSettings/projectESL_output.nc')
gesla3_99 = xr.open_dataset('/Users/timhermans/Documents/GitHub/projectESL/output/COAST_RP/GESLA3_99pct/projectESL_output.nc')


locs_intersect = np.intersect1d(gtsm_ref.locations,crp.locations)
crp = crp.sel(locations = locs_intersect)
gtsm_ref = gtsm_ref.sel(locations = locs_intersect)
gtsm_99 = gtsm_99.sel(locations = locs_intersect)
'''
gtsm_ref = gtsm_ref.sel(locations = locs_intersect)
gtsm_99 = gtsm_99.sel(locations = locs_intersect)

gesla3_ref = gesla3_ref.sel(locations = locs_intersect)
gesla3_99 = gesla3_99.sel(locations = locs_intersect)
'''
'''
fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(b.lon,b.lat,
                c=b.z_hist.sel(f=1).isel(qnt=1)-c.z_hist.sel(f=1).isel(qnt=1),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-.5,vmax=.5)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('Muis23 without minus Muis2023 with preprocessing')
'''


fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(gesla3_ref.lon,gesla3_ref.lat,
                c=100*((gesla3_99.z_hist.sel(f=.01).isel(qnt=1))-(gesla3_ref.z_hist.sel(f=.01).isel(qnt=1)))/(gesla3_ref.z_hist.sel(f=.01).isel(qnt=1)),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-10,vmax=10)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [%]')
ax.set_title('gesla3 ref minus 99pct')


fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(gtsm_ref.lon,gtsm_ref.lat,
                c=100*((gtsm_99.z_hist.sel(f=.01).isel(qnt=1))-(gtsm_ref.z_hist.sel(f=.01).isel(qnt=1)))/(gtsm_ref.z_hist.sel(f=.01).isel(qnt=1)),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-10,vmax=10)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=100yr) [%]')
ax.set_title('gtsm ref minus 99pct')






fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(gesla3_ref.lon,gesla3_ref.lat,
                c=100*((gesla3_99.z_hist.sel(f=.01).isel(qnt=1)-gesla3_99.z_hist.sel(f=1).isel(qnt=1))-(gesla3_ref.z_hist.sel(f=.01).isel(qnt=1)-gesla3_ref.z_hist.sel(f=1).isel(qnt=1)))/(gesla3_ref.z_hist.sel(f=.01).isel(qnt=1)-gesla3_ref.z_hist.sel(f=1).isel(qnt=1)),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-10,vmax=10)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('gesla3 ref minus 99pct')


fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(gtsm_ref.lon,gtsm_ref.lat,
                c=100*((gtsm_99.z_hist.sel(f=.01).isel(qnt=1)-gtsm_99.z_hist.sel(f=1).isel(qnt=1))-(gtsm_ref.z_hist.sel(f=.01).isel(qnt=1)-gtsm_ref.z_hist.sel(f=1).isel(qnt=1)))/(gtsm_ref.z_hist.sel(f=.01).isel(qnt=1)-gtsm_ref.z_hist.sel(f=1).isel(qnt=1)),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-20,vmax=20)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('gtsm ref minus 99pct')



fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(gtsm_ref.lon,gtsm_ref.lat,
                c=100*((crp.z_hist.sel(f=.01).isel(qnt=1))-(gtsm_99.z_hist.sel(f=.01).isel(qnt=1)))/(gtsm_99.z_hist.sel(f=.01).isel(qnt=1)),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-20,vmax=20)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=100yr) [m]')
ax.set_title('coast-rp minus gtsm')



'''

fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(a.lon,a.lat,
                c=a.z_hist.sel(f=1).isel(qnt=1),
                transform=ccrs.PlateCarree(),s=3,cmap='YlGnBu',vmin=0,vmax=3)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('Muis23 without minus Muis2023 with preprocessing')



fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(a0.lon,a0.lat,
                c=a0.z_hist.sel(f=1).isel(qnt=1),
                transform=ccrs.PlateCarree(),s=3,cmap='YlGnBu',vmin=0,vmax=3)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('Muis23 without minus Muis2023 with preprocessing')



fig=plt.figure(figsize=(9,9.5)) #generate figure 
gs = fig.add_gridspec(1,1)
gs.update(hspace=.2)

ax = plt.subplot(gs[0,0],projection=ccrs.Robinson(central_longitude=0))

ax.add_feature(cartopy.feature.OCEAN, zorder=0,facecolor='grey')
ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor='grey')

sc = ax.scatter(a0.lon,a0.lat,
                c=a0.z_hist.sel(f=.01).isel(qnt=1)-a.z_hist.sel(f=.01).isel(qnt=1),
                transform=ccrs.PlateCarree(),s=3,cmap='coolwarm',vmin=-.1,vmax=.1)
cax=ax.inset_axes(bounds=(0, -.1,1,.075))
cb=fig.colorbar(sc, cax=cax,orientation='horizontal',label='Difference in return height (RP=1yr) [m]')
ax.set_title('GTSM minus Vousdoukas2018')
'''