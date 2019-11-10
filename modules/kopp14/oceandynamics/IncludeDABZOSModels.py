import numpy as np
import os
import re
import glob
import scipy.io
import scipy.sparse
import h5py
import time

''' IncludeDABZOSModels.py

This script parses through a directory of models to determine if the output from each
model should be included in the analysis.

Parameters: 
model_dir = Directory of model output. Each model is a subdirectory within this one.
varnames =  Names of the variables of interest. Will be subdirectories of individual 
			model directories.
years = Years of interest.
site_ids = Site ID numbers to extract from the ZOS file.

Return: 
model_list = Vector of model names that are to be included (nmodels)
ZOS = Sea-level change due to ocean dynamics from CMIP5 output (years, nmodels, sites)

'''

def angd(lat0, lon0, lat, lon):
	
	# Convert the input from degrees to radians
	(lat0, lon0) = np.radians((lat0, lon0))
	(lat, lon) = np.radians((lat, lon))
	
	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + \
	(np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),\
	(np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))
	
	# Convert the results from radians to degrees and return
	return(np.degrees(temp))

def CalcWeights(qlat, qlon, lat, lon, idwrad, idwpow, idwmin):
	
	# Calculate the distance between the query lat/lon and the grid points
	ad = angd(lat, lon, qlat, qlon)
	ad[np.nonzero(ad <= idwmin)] = idwmin
	
	# Find the indices that are below the IDW threshold
	idw_idx = np.nonzero(ad <= idwrad)
	
	# Calculate the weights
	idw_weights = ad[idw_idx]**-idwpow
	
	return(idw_idx, idw_weights)
	
def IDW(val, idw_weights, idw_idx):

	# Apply the weights to the values
	qval = np.nansum(val[idw_idx[0], idw_idx[1],:] * idw_weights[:,np.newaxis], axis=0) / np.nansum(idw_weights)
	
	return(qval)

def IncludeDABZOSModels(model_dir, rcp_scen, focus_site_lats, focus_site_lons):
	
	# Initialize the model list and data arrays
	model_list = []
	ZOS = []
	
	# Load the model lat/lon files and extract the data
	mat_file = os.path.join(model_dir, "slrlat.mat")
	mat = scipy.io.loadmat(mat_file)
	lats = mat['nlat']
	
	mat_file = os.path.join(model_dir, "slrlon.mat")
	mat = scipy.io.loadmat(mat_file)
	lons = mat['nlon']
	
	#start_time = time.clock()
	
	# Initialize variables to hold the IDW weights and indices
	all_idx = []
	all_weights = []
	
	# Initialize IDW parameters
	idw_rad = 3.5
	idw_pow = 3.0
	idw_min = 0.005
	
	# Calculate the weights for all the sites
	# Note: this assumes grid for each model is the same
	for i in np.arange(len(focus_site_lats)):
		
		# Calculate the weights
		(idx, weights) = CalcWeights(focus_site_lats[i], focus_site_lons[i], lats, lons, idw_rad, idw_pow, idw_min)
		all_idx.append(idx)
		all_weights.append(weights)
	
	# Loop through available models in the model/rcp directory
	rcp_dir = os.path.join(model_dir, rcp_scen)
	for model in glob.iglob(os.path.join(rcp_dir, "*zos.mat")):
		
		# Get the model name from the file name
		match_string = rcp_dir + r"/(.*)_rcp\d\d_zos.mat$"
		model_match = re.search(match_string, model)
		model_name = model_match.group(1).lower()
		
		# Append this model to the model list
		model_list.append(model_name)
		
		# Load model file
		with h5py.File(model, 'r') as mat:
			raw_zos = np.array(mat['zosv2']).T
		
		# Calculate the zos values for all sites from this model
		model_zos = map(lambda idx, w: IDW(raw_zos, w, idx), all_idx, all_weights)
		
		# Remove the 1860 value to be consistent with the ZOSTOGA data and pad the year 
		# dimension to year 2300 to match length of ZOSTOGA (440 - 1 = 439 entries)
		model_zos = map(lambda this_zos: (list(this_zos) + 440 * [np.nan])[1:440], model_zos)
		ZOS.append(model_zos)
		
	# Convert the ZOS list to a numpy array
	ZOS = np.array(ZOS)
	
	# Permute the ZOS array (years, models, sites)
	ZOS = np.transpose(ZOS, axes=(2,0,1))
		
	return(model_list, ZOS)

if __name__ == "__main__":
	
	start_time = time.clock()
	site_lats = [40.1, 32.6, 40.1, 40.1]
	site_lons = [-72.3, -117.8, -72.3, -72.3]
	model_dir = "/Users/ggg46/research/slr_framework/data/kopp14/zosfinal"
	rcp_scen = "rcp85"
	(modellist, zos) = IncludeDABZOSModels(model_dir, rcp_scen, site_lats, site_lons)
	end_time = time.clock()
	
	print(end_time - start_time)
	print(zos[238,0,:])