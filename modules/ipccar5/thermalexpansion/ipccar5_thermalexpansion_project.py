# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19

import os
import numpy as np
import pickle
import argparse
import time
from netCDF4 import Dataset


class ProjectionError(Exception):
	pass
	

def ar5_project_thermalexpansion(rng_seed, pyear_start, pyear_end, pyear_step, nsamps, pipeline_id):
	
	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	
	# Load the preprocessed data
	data_file = "{}_data.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")
	
	# Extract the data variables
	my_data = pickle.load(f)
	f.close()
	
	exp_mean = my_data['exp_mean']
	exp_sd = my_data['exp_sd']
	data_years = my_data['data_years']
	startyr = my_data["startyr"]
	scenario = my_data['scenario']
	
	# Subset the data to the target years
	year_idx = np.isin(data_years, targyears)
	exp_mean = exp_mean[year_idx]
	exp_sd = exp_sd[year_idx]
	data_years = data_years[year_idx]
	
	# Set the seed for the random number generator
	rng = np.random.default_rng(rng_seed)

	# Generate perfectly correlated samples
	z=rng.standard_normal(nsamps)[:,np.newaxis]
	
	# For each quantity, mean + standard deviation * normal random number
	zx=exp_mean + (exp_sd * z)
	
	# Number of years in the data record
	nyr = len(data_years)
	
	# Save the global thermal expansion projections to a pickle
	output = {"zx": zx, "data_years": data_years, "startyr": startyr, "scenario": scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Convert to mm
	zx = zx * 1000  # Convert to mm

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(data_years))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim  = rootgrp.createDimension("locations", 1)
	
	# Populate dimension variables
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))
	loc_var  = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var  = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var  = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from thermal expansion according to AR5 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, scenario)
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:]  = targyears
	samp_var[:]  = np.arange(nsamps)
	samps[:,:,:] = zx[:,:,np.newaxis]
	lat_var[:] 	 = np.inf
	lon_var[:] 	 = np.inf
	loc_var[:] 	 = -1


	# Close the netcdf
	rootgrp.close()	

	return(0)

	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glacier projection stage for the AR5 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate [default=1000]", default=1000, type=int)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	ar5_project_thermalexpansion(args.seed, args.pyear_start, args.pyear_end, args.pyear_step, args.nsamps, args.pipeline_id)
	
	exit()
