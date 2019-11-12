import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from netCDF4 import Dataset
from read_bkgdrate import read_bkgdrate
from AssignFP import AssignFP


''' kopp14_postprocess_glaciers.py

This script runs the glacier post-processing task for the Kopp 2014 workflow. This task
uses the global projections from the 'kopp14_project_glaciers' script and applies
spatially resolved fingerprints to the glacier and ice cap contributions. The result is a 
netCDF4 file that contains spatially and temporally resolved samples of GIC contributions 
to local sea-level rise

Parameters: 
focus_site_ids = Location IDs for localization (from PSMSL)
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing local contributions from GIC

'''

def kopp14_postprocess_glaciers(focus_site_ids, pipeline_id):
	
	# Read in the global projections
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)
	
	# Extract the projection data from the file
	my_data = pickle.load(f)
	gicsamps = my_data["gicsamps"]
	f.close()
	
	# Read in the fingerprint information
	fpfile = "{}_fp.pkl".format(pipeline_id)
	try:
		f = open(fpfile, 'rb')
	except:
		print("Cannot open fpfile\n")
		sys.exit(1)
	
	# Extract the fingerprint information from the file
	my_data = pickle.load(f)
	fpmapperids = my_data["fpmapperids"]
	fpmaps = my_data["fpmaps"]
	f.close()
	
	# Read in the configuration information
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configfile\n")
		sys.exit(1)
	
	# Extract the fingerprint information from the file
	my_data = pickle.load(f)
	targyears = my_data["targyears"]
	rcp_scenario = my_data["rcp_scenario"]
	f.close()
	
	
	# Load the site locations	
	ratefile = os.path.join(os.path.dirname(__file__), "bkgdrate.tsv")
	(_, site_ids, site_lats, site_lons) = read_bkgdrate(ratefile, True)
	
	# FOR SIMPLICITY, LOCALIZE TO ONLY A FEW LOCATIONS
	if np.any([x >= 0 for x in focus_site_ids]):
		_, _, site_inds = np.intersect1d(focus_site_ids, site_ids, return_indices=True)
		site_ids = site_ids[site_inds]
		site_lats = site_lats[site_inds]
		site_lons = site_lons[site_inds]
	
	# Initialize variable to hold the localized projections
	(nsamps, nregions, ntimes) = gicsamps.shape
	nsites = len(site_ids)
	local_sl = np.full((nsites, nsamps, ntimes), 0.0) 
	
	# Loop through the GIC regions
	for i in np.arange(0,nregions):
		
		# Get the fingerprint file name for this region
		thisRegion = fpmaps[i]
		
		# Get the fingerprints for these sites from this region
		regionfile = os.path.join(os.path.dirname(__file__), "FPRINT", "fprint_{0}.mn".format(thisRegion))
		regionfp = AssignFP(regionfile, site_lats, site_lons)
		
		# Multiply the fingerprints and the projections and add them to the running total
		# over the regions
		local_sl = local_sl + np.transpose(np.multiply.outer(gicsamps[:,i,:], regionfp), (2,0,1))
	
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
	rootgrp.description = "Local SLR contributions from glaciers and ice caps according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: Kopp 2014 GIC workflow - {0}".format(rcp_scenario)
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
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the GIC component of the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--site_ids', help="Site ID numbers (from PSMSL database) to make projections for")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
		
	# Parse the arguments
	args = parser.parse_args()
	
	# Convert the string of site_ids to a list
	site_ids = [int(x) for x in re.split(",\s*", str(args.site_ids))]
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14_postprocess_glaciers(site_ids, args.pipeline_id)
	
	# Done
	exit()