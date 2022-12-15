import numpy as np
import pickle
import sys
import os
import argparse
import time
import re
from read_locationfile import ReadLocationFile
from netCDF4 import Dataset

''' ar5_postprocess_thermalexp.py

This script runs the thermal expansion postprocessing task for the AR5 workflow. 
This task generates localized contributions to sea-level change due to thermosteric heating
of the ocean.

Parameters: 
locationfile = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing the local sea-level rise projections

Note: The AR5 workflow does not apply a fingerprint to the thermal expansion component. 
As such, this script applies an implicit fingerprint coefficient of '1' and produces
localized slr projections for the site ids provided.

'''

def ar5_postprocess_thermalexp(locationfilename, pipeline_id):
	
	# Load the projection file
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projection file {}\n").format(projfile)
	
	# Extract the configuration variables
	my_proj = pickle.load(f)
	f.close()
	
	targyears = my_proj["data_years"]
	thermsamps = my_proj["zx"]
	startyr = my_proj["startyr"]
	rcp_scenario = my_proj['scenario']
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
	# Initialize variable to hold the localized projections
	(nsamps, ntimes) = thermsamps.shape
	nsites = len(site_ids)
	
	# Apply the effective fingerprint of "1" to the global projections for each site
	local_sl = np.tile(thermsamps, (nsites, 1, 1))
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", ntimes)
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim  = rootgrp.createDimension("locations", nsites)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var  = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var  = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var  = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from thermal expansion according to AR5 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}, Start year = {2}".format(pipeline_id, rcp_scenario, startyr)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:]  = targyears
	samp_var[:]  = np.arange(nsamps)
	samps[:,:,:] = np.transpose(local_sl[:,:,:],(1,2,0))
	lat_var[:] 	 = site_lats
	lon_var[:] 	 = site_lons
	loc_var[:] 	 = site_ids

	# Close the netcdf
	rootgrp.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion postprocessing stage for the AR5 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the AR5 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	ar5_postprocess_thermalexp(args.locationfile, args.pipeline_id)
	
	exit()