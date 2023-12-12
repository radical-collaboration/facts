import numpy as np
import argparse
import pickle
import sys
import os
from netCDF4 import Dataset
import time

''' dp16_project_icesheet.py

Runs the dp16 icesheet projection stage.

Parameters: 
nsamps              Number of samples to produce
replace             Allow sampling with replacement
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def dp16_project_icesheet(nsamps, pyear_start, pyear_end, pyear_step, pipeline_id, replace, rngseed):
	
	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)
	
	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["years"]
	wais = my_data["wais_samps"]
	eais = my_data["eais_samps"]
	scenario = my_data["scenario"]
	
	# Define the target projection years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	
	# Extract the pool size from the data
	pool_size = eais.shape[1]
	
	# Find the data years that overlap with the target projection years
	(_, datayr_idx, targyear_idx) = np.intersect1d(years, targyears, return_indices = True)
	
	# Generate the sample indices
	rng = np.random.default_rng(rngseed)
	sample_idx = rng.choice(pool_size, size=nsamps, replace=replace)
	
	# Store the samples
	wais_samps = wais[datayr_idx[:,np.newaxis],sample_idx[np.newaxis,:]]
	eais_samps = eais[datayr_idx[:,np.newaxis],sample_idx[np.newaxis,:]]
	ais_samps = wais_samps + eais_samps
	   
    # Store the variables in a pickle
	output = {'eais_samps': eais_samps, 'wais_samps': wais_samps, \
				'scenario': scenario, 'targyears': targyears[targyear_idx]}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Write the global projections to output netCDF files
	WriteNetCDF(eais_samps, "EAIS", targyears[targyear_idx], scenario, pipeline_id)
	WriteNetCDF(wais_samps, "WAIS", targyears[targyear_idx], scenario, pipeline_id)
	WriteNetCDF(ais_samps, "AIS", targyears[targyear_idx], scenario, pipeline_id)
	
	return(None)


def WriteNetCDF(icesamps, icetype, data_years, scenario, pipeline_id):

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_{1}_globalsl.nc".format(pipeline_id, icetype))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	nyr = len(data_years)
	nsamps = icesamps.shape[1]
	year_dim = rootgrp.createDimension("years", nyr)
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from {} according to DP16 workflow".format(icetype)
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, scenario)
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:]  = data_years
	samp_var[:]  = np.arange(nsamps)
	samps[:,:,:] = np.transpose(icesamps[:,:,np.newaxis],(1,0,2))
	lat_var[:] 	 = np.inf
	lon_var[:] 	 = np.inf
	loc_var[:] 	 = -1

	# Close the netcdf
	rootgrp.close()	

	return(0)


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the generic module \'directsample\' projection stage.",\
	epilog="Note: This is meant to be run as part of the \'genmod\' module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 10000)", default=10000, type=int)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--replace', help="Allow sampling with replacement [default = 1]", choices=('0','1'), default=1)
	parser.add_argument('--seed', help="Seed for the random number generator [default = 1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection stage with the provided arguments
	dp16_project_icesheet(args.nsamps, args.pyear_start, args.pyear_end, args.pyear_step, args.pipeline_id, args.replace, args.seed)
	
	exit()