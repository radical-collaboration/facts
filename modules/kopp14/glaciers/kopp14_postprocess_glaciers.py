import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from netCDF4 import Dataset
from read_locationfile import ReadLocationFile
from AssignFP import AssignFP


''' kopp14_postprocess_glaciers.py

This script runs the glacier post-processing task for the Kopp 2014 workflow. This task
uses the global projections from the 'kopp14_project_glaciers' script and applies
spatially resolved fingerprints to the glacier and ice cap contributions. The result is a 
netCDF4 file that contains spatially and temporally resolved samples of GIC contributions 
to local sea-level rise

Parameters: 
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing local contributions from GIC

'''

def kopp14_postprocess_glaciers(locationfilename, pipeline_id):
	
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
	baseyear = my_data["baseyear"]
	f.close()
	
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
	# Initialize variable to hold the localized projections
	(nsamps, nregions, ntimes) = gicsamps.shape
	nsites = len(site_ids)
	local_sl = np.full((nsites, nsamps, ntimes), 0.0) 
	
	# Loop through the GIC regions
	for i in np.arange(0,nregions):
		
		# Get the fingerprint file name for this region
		thisRegion = fpmaps[i]
		
		# Get the fingerprints for these sites from this region
		regionfile = os.path.join(os.path.dirname(__file__), "FPRINT", "fprint_{0}.nc".format(thisRegion))
		regionfp = AssignFP(regionfile, site_lats, site_lons)
		
		# Multiply the fingerprints and the projections and add them to the running total
		# over the regions
		local_sl = local_sl + np.transpose(np.multiply.outer(gicsamps[:,i,:], regionfp), (2,0,1))
	
	# Calculate the quantiles
	out_q = np.unique(np.append(np.linspace(0,1,101), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))
	nq = len(out_q)
	local_sl_q = np.nanquantile(local_sl, out_q, axis=1)
		
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
	localslq = rootgrp.createVariable("localSL_quantiles", "i2", ("quantiles", "nsites", "years"), zlib=True, complevel=4)
	#localslq.scale_factor = 0.1

	# Assign attributes
	rootgrp.description = "Local SLR contributions from glaciers and ice caps according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: Kopp 2014 GIC workflow - {0}; Base year {1}".format(rcp_scenario, baseyear)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	localslq.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	q_var[:] = out_q
	localslq[:,:,:] = local_sl_q

	# Close the netcdf
	rootgrp.close()
	
if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the GIC component of the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
		
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14_postprocess_glaciers(args.locationfile, args.pipeline_id)
	
	# Done
	exit()
