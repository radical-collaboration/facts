import numpy as np
import pickle
import os
import time
import argparse
from netCDF4 import Dataset

''' gilford_postprocess_icesheets.py

This script runs the icesheet post processing task for the Gilford icesheet emulator. 

Parameters: 
config_file = Configuration file produced in the preprocessing task
proj_file = Pickled file with configuration information generated in pre-processing

Output: NetCDF file containing projections

Note:

'''

def gilford_postprocess_icesheets(proj_file):
	
	# Load the projection file
	try:
		f = open(proj_file, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_proj = pickle.load(f)
	f.close()
	
	proj_years = my_proj["years"]
	proj_samps = my_proj["samps"] * 1000   # Convert to mm
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gilford_output.nc"), "w", format="NETCDF4")

	# Define Dimensions
	#site_dim = rootgrp.createDimension("nsites", len(sites))
	year_dim = rootgrp.createDimension("years", len(proj_years))
	samp_dim = rootgrp.createDimension("samples", proj_samps.shape[1])

	# Populate dimension variables
	#lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	#lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	#id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	ais_proj = rootgrp.createVariable("AIS", "f4", ("years", "samples"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Sea-level rise contributions from AIS according to Gilford Emulator"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "SLR Framework: Gilford Emulator"
	#lat_var.units = "Degrees North"
	#lon_var.units = "Degrees West"
	#id_var.units = "[-]"
	year_var.units = "[-]"
	samp_var.units = "[-]"
	ais_proj.units = "mm"

	# Put the data into the netcdf variables
	#lat_var[:] = site_lats
	#lon_var[:] = site_lons
	#id_var[:] = site_ids
	year_var[:] = proj_years
	samp_var[:] = np.arange(0,proj_samps.shape[1])
	ais_proj[:,:] = proj_samps

	# Close the netcdf
	rootgrp.close()


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the ice sheet post-processing stage for the Gilford Ice Sheet Emulator workflow",\
	epilog="Note: This is meant to be run as part of the Gilford Ice Sheet module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--proj_file', help="Projection file generated from the projection stage",\
	default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gilford_icesheets_projections.pkl"))
	
	# Parse the arguments
	args = parser.parse_args()
			
	# Run the post-processing stage
	gilford_postprocess_icesheets(args.proj_file)
	
	exit()