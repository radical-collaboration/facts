import numpy as np
import pickle
import sys
import os
import argparse
import time
import re
from read_locationfile import ReadLocationFile
from netCDF4 import Dataset

''' kopp14_postprocess_thermalexp.py

This script runs the thermal expansion postprocessing task for the Kopp 2014 workflow. 
This task generates localized contributions to sea-level change due to thermosteric heating
of the ocean.

Parameters: 
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing the local sea-level rise projections

Note: The Kopp14 workflow does not apply a fingerprint to the thermal expansion component. 
As such, this script applies an implicit fingerprint coefficient of '1' and produces
localized slr projections for the site ids provided.

'''

def kopp14_postprocess_thermalexp(locationfilename, pipeline_id):
	
	# Load the configuration file
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projection file {}\n").format(projfile)
	
	# Extract the configuration variables
	my_proj = pickle.load(f)
	f.close()
	
	targyears = my_proj["targyears"]
	thermsamps = my_proj["thermsamps"]
	rcp_scenario = my_proj['rcp_scenario']
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
	# Initialize variable to hold the localized projections
	(nsamps, ntimes) = thermsamps.shape
	nsites = len(site_ids)
	
	# Apply the effective fingerprint of "1" to the global projections for each site
	local_sl = np.tile(thermsamps, (nsites, 1, 1))
	
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
	rootgrp.description = "Local SLR contributions from thermal expansion according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, rcp_scenario)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
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
	parser = argparse.ArgumentParser(description="Run the thermal expansion postprocessing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	kopp14_postprocess_thermalexp(args.locationfile, args.pipeline_id)
	
	exit()