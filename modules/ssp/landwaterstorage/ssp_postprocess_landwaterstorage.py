import numpy as np
import pickle
import sys
import os
import argparse
import time
import re
from read_locationfile import ReadLocationFile
from AssignFP import AssignFP
from netCDF4 import Dataset

''' ssp_postprocess_landwaterstorage.py

This script runs the land water storage postprocessing task from the SSP module set. 
This task generates localized contributions to sea-level change due to land water storage.

Parameters: 
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing the local sea-level rise projections

'''

def ssp_postprocess_landwaterstorage(locationfilename, pipeline_id):
	
	# Load the configuration file
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projection file {}\n").format(projfile)
	
	# Extract the configuration variables
	my_proj = pickle.load(f)
	f.close()
	
	targyears = my_proj["years"]
	scen = my_proj['scen']
	lwssamps = np.transpose(my_proj["lwssamps"])
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)	
	
	# Initialize variable to hold the localized projections
	(nsamps, ntimes) = lwssamps.shape
	nsites = len(site_ids)
	
	# Apply the fingerprints
	fpfile = os.path.join(os.path.dirname(__file__), "REL_GROUNDWATER_NOMASK.nc")
	fpsites = AssignFP(fpfile, site_lats, site_lons)
	local_sl = lwssamps[np.newaxis,:,:] * fpsites[:,np.newaxis,np.newaxis]
	
	# Calculate the quantiles
	out_q = np.unique(np.append(np.linspace(0,1,101), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))
	nq = len(out_q)
	local_sl_q = np.nanquantile(local_sl, out_q, axis=1)
	
	# Calculate the mean and sd of the samples
	local_sl_mean = np.nanmean(local_sl, axis=1)
	local_sl_sd = np.nanstd(local_sl, axis=1)
		
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", nsites)
	year_dim = rootgrp.createDimension("years", ntimes)
	q_dim = rootgrp.createDimension("quantiles", nq)

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable("localSL_quantiles", "f4", ("quantiles", "nsites", "years"), zlib=True, least_significant_digit=2)
	localslmean = rootgrp.createVariable("localSL_mean", "f4", ("nsites", "years"), zlib=True, least_significant_digit=2)
	localslsd = rootgrp.createVariable("localSL_std", "f4", ("nsites", "years"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from land water storage from the SSP module set"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, scen)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	localslq.units = "mm"
	localslmean.units = "mm"
	localslsd.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	q_var[:] = out_q
	localslq[:,:,:] = local_sl_q
	localslmean[:,:] = local_sl_mean
	localslsd[:,:] = local_sl_sd

	# Close the netcdf
	rootgrp.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the land water storage postprocessing stage from the SSP module set",\
	epilog="Note: This is meant to be run as part of the SSP module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	ssp_postprocess_landwaterstorage(args.locationfile, args.pipeline_id)
	
	exit()