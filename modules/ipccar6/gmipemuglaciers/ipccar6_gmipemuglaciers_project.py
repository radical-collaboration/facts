import numpy as np
import argparse
import pickle
import os
import re
import time
import sys
from netCDF4 import Dataset
from scipy.stats import norm

''' ipccar6_project_gmipemuglaciers.py

This is the projection stage for the GMIP2 emulated glaciers component of the IPCC AR6 module set.

Parameters:
nsamps              Number of samples to produce
replace             Allow sampling with replacement
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline


'''

def ipccar6_project_gmipemuglaciers(nsamps, pipeline_id, replace, rngseed):

	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["targyears"]
	scenario = my_data["scenario"]
	gic_samples = my_data["gic_samps"]
	baseyear = my_data["baseyear"]
	model_driver = my_data["model_driver"]

	# Generate the sample indices
	rng = np.random.default_rng(rngseed)
	sample_inds = rng.choice(gic_samples.shape[1], size=nsamps, replace=replace)

	# Extract the samples
	gic_region_samps = gic_samples[:,sample_inds,:]

	# Sum up the regions
	gic_samps = np.sum(gic_region_samps, axis=0)

    # Store the variables in a pickle
	output = {'gic_samps': gic_region_samps, 'years': years, 'scenario': scenario, \
				'baseyear': baseyear, 'model_driver': model_driver}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the projections to the netCDF files
	WriteNetCDF(pipeline_id, gic_samps, years, nsamps, scenario, baseyear)

	return(0)


def WriteNetCDF(pipeline_id, global_samps, years, nsamps, scenario, baseyear):

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(years))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = "Global SLR contribution from glaciers from the GMIP2 emulated workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}. ".format(pipeline_id)
	rootgrp.baseyear = baseyear
	rootgrp.scenario = scenario
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = years
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = global_samps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()

if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the GMIP2 emulated glaciers projection stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 10)", default=10, type=int)
	parser.add_argument('--replace', help="Allow sampling with replacement (default = 1)", choices=(0,1), type=int, default=1)
	parser.add_argument('--seed', help="Seed for the random number generator (default = 1234)", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the projection stage with the provided arguments
	ipccar6_project_gmipemuglaciers(args.nsamps, args.pipeline_id, args.replace, args.seed)

	exit()
