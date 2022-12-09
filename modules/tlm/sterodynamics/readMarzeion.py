import numpy as np
import os
import glob
import re

''' readMarzeion.py

This script loads the glacier volumetric change data

Parameters: 
scen = RCP scenario (e.g. "rcp85")
glacdir = Full path to the directory that contains the glacier data
fpmap = Full path and file name that contains the fingerprint region information
discardAntarctica = Whether or not to discard Antarctica from the returned values (default: False)

Return: 
GIC = Array of glacier contributions to sea-level rise (year, region, model)
GICse = Standard error for each entry in GIC (year, region, model)
GICyrs = Years associated with the glacier data (year, model)
GICmodel = Names of models from which data were collected (model)
fpmapperids = Fingerprint Mapper IDs (IceID from Marzeion file)
fpmaps = Fingerprint IDs (that match fingerprint data files)
GICnames = Names of the glacier regions

Note: 

'''

def readMarzeion(scen, glacdir, fpmap, discardAntarctica=False):
	
	# Set the baseyear of interest
	baseyear = 2005
	
	# Open the file
	fpmap_data = np.genfromtxt(fpmap, dtype=None, names=True, delimiter=',', encoding=None)
	
	# Extract the data
	fpmapperids = fpmap_data['IceID']
	fpmaps = fpmap_data['FPID']
	GICnames = fpmap_data['IceName']
	
	# Discard Antarctica if necessary
	region_inds = np.arange(0,len(fpmapperids))
	if(discardAntarctica):
		region_inds = np.nonzero(GICnames != "AntarcticGIC")[0]
		fpmapperids = fpmapperids[region_inds]
		fpmaps = fpmaps[region_inds]
		GICnames = GICnames[region_inds]
	
	# Get a list of regional accumulation files to parse through
	glob_string = os.path.join(glacdir, "dV_regional_accum*" + scen + "*.txt")
	files = glob.glob(glob_string)
	
	# Initialize variables to hold glacier data
	# Note: It appears each file may contain data from year 1850 - 2300
	dim_year = 2300 - 1850 + 1
	dim_region = len(region_inds)
	dim_model = len(files)
	GIC = np.full((dim_year, dim_region, dim_model), np.nan)
	GICse = np.full((dim_year, dim_region, dim_model), np.nan)
	GICyrs = np.full((dim_year, dim_model), np.nan)
	GICmodel = []
	
	# Loop through the files
	for i in np.arange(0,len(files)):
		
		# Import the data from the file
		glac_data = np.genfromtxt(files[i], dtype="float")
		
		# Find the row that represents the baseyear
		baseyear_ind = np.nonzero(glac_data[:,0] == baseyear)[0]
		
		# Center the data on the baseyear
		glac_data[:,1:] = glac_data[:,1:] - glac_data[baseyear_ind,1:]
		
		# Extract the data
		GIC[:,:,i] = glac_data[:,region_inds+1]
		GICse[:,:,i] = glac_data[:,region_inds+20]
		GICyrs[:,i] = glac_data[:,0]
		
		# Extract the name of the model from the file name
		regex_test = re.search(r"accum_(.+)_rcp", files[i])
		GICmodel.append(regex_test.group(1).lower())
	
	# Return the objects
	return(GIC, GICse, GICyrs, GICmodel, fpmapperids, fpmaps, GICnames)