'''
Utility functions for projectESL. 

@author: Tim Hermans
t(dot)h(dot)j(dot)hermans@uu(dot)nl
'''

import numpy as np
import xarray as xr
import os

def if_scalar_to_list(variable):
    if np.isscalar(variable):
        variable_as_list = [variable]
        return variable_as_list
    else:
        return variable

def kmdist(qlat,qlon,lats,lons):
    qlat = np.array(qlat)
    qlon = np.array(qlon)
    lats = np.array(lats)
    lons = np.array(lons)
    distances = 2*np.arcsin( np.sqrt(
            np.sin( (np.pi/180) * 0.5*(lats-qlat) )**2 +
            np.cos((np.pi/180)*qlat)*np.cos((np.pi/180)*lats)*np.sin((np.pi/180)*0.5*(lons-qlon))**2) )
    
    return distances*6371

def angdist(lat0, lon0, lat, lon):
	'''calculate angular distance between coordinates (lat0,lon0) and (lat,lon)'''
	# Convert the input from degrees to radians
	(lat0, lon0) = np.radians((lat0, lon0))
	(lat, lon) = np.radians((lat, lon))
	
	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + \
	(np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),\
	(np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))
	
	# Convert the results from radians to degrees and return
	return(np.degrees(temp))

def mindist(qlat, qlon, lats, lons, limit=10):
	'''find indices of the smallest distances between (qlat,qlon) and (lats,lons) within kilometric distance limit "limit" '''
	# Calculate the angular distance
	dist = kmdist(lats, lons, qlat, qlon)
	
	# If the minimum distance is beyond the limit, print a warning and return None
	if(np.amin(dist) > limit):
		return(None)
	
	else:
		# Perform an indirect sort of the distances
		sort_idx = np.argsort(dist)
		
		# Find which ones fall within the radius limit
		min_idx = sort_idx[np.flatnonzero(dist[sort_idx] <= limit)]
		
		return(min_idx)
    
def download_ar6_full_sample_projections_at_tgs(wf,ssp,out_dir): #~1.2Gb!
    #see https://github.com/Rutgers-ESSP/IPCC-AR6-Sea-Level-Projections
    ds_url = ('https://storage.googleapis.com/ar6-lsl-simulations-public-standard/'
         'tide-gauges/full_sample_workflows/'+wf+'/'+ssp+'/total-workflow.zarr'
       )
    ds = xr.open_dataset(ds_url, engine='zarr', chunks='auto')
    return ds.to_zarr(os.path.join(out_dir,'full_sample_total_'+wf+'_'+ssp+'.zarr'),mode='w')

def download_ar6_full_sample_projections_gridded(wf,ssp,out_dir): #~52Gb!
    #see https://github.com/Rutgers-ESSP/IPCC-AR6-Sea-Level-Projections
    ds_url = ('https://storage.googleapis.com/ar6-lsl-simulations-public-standard/'
         'gridded/full_sample_workflows/'+wf+'/'+ssp+'/total-workflow.zarr'
       )
    ds = xr.open_dataset(ds_url, engine='zarr', chunks='auto')
    return ds.to_zarr(os.path.join(out_dir,'full_sample_total_'+wf+'_'+ssp+'.zarr'),mode='w')

def add_ar6_full_sample_projections_to_locations(input_locations,slr_fn,nsamps,period):
    '''adds nearest sea-level projection samples from AR6 (in file 'slr_fn') to dataset with lon/lat coordinates at locations'''
    ds = xr.open_dataset(slr_fn, engine='zarr', chunks='auto') #open projections
    ds_stacked = ds.stack(locations=['lon','lat']) #stack along lon/lat
    
    slr_lons = ds_stacked.lon.values
    slr_lats = ds_stacked.lat.values
    
    mask = np.isnan(ds.isel(years=0,samples=0).sea_level_change).stack(locations=['lon','lat']).values #get land mask
    min_idx = [np.argmin(kmdist(qlat,qlon,slr_lats,slr_lons)+999*mask) for qlat,qlon in zip(input_locations.lat.values, input_locations.lon.values)] #get nearest projections to input_locations, don't use land

    ds_at_input_locations = ds_stacked.isel(locations=min_idx).drop('lat') #drop multiindex
    ds_at_input_locations['locations']=input_locations['locations'].values #assign locations from input_locations file
    
    input_locations['sea_level_change'] = ds_at_input_locations['sea_level_change'] #add SLC variable to input_locations file
    input_locations = input_locations.sel(years=slice(str(period[0]),str(period[1]))).isel(samples=np.arange(nsamps)) #select requested number of samples (samples are randomly ordered in the input) and period

    return input_locations
