import numpy as np
import pickle
import sys
import os
import argparse
import time
from netCDF4 import Dataset
from scipy.stats import norm
from scipy.stats import t

''' kopp14_project_glaciers.py

This script runs the glacier projection task for the Kopp 2014 workflow. 
This task generates global contributions to sea-level change due to glacier and ice cap
melt.

Parameters: 
pipeline_id = Unique identifier for the pipeline running this code
nsamps = Numer of samples to produce
seed = Seed for the random number generator

Output: 
- Pickle file containing Glacier and Ice Cap contributions to sea-level rise (nsamps, regions, years)
- NetCDF files containing the sum total global contributions across all glacier regions

'''

def kopp14_project_glaciers(nsamps, seed, pipeline_id):
	
	# Load the fit file
	fitfile = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(fitfile, 'rb')
	except:
		print("Cannot open fit file\n")
	
	# Extract the fit variables
	my_fit = pickle.load(f)
	f.close()
	
	meanGIC = my_fit["meanGIC"]
	T = my_fit["T"]	
	NGIC = my_fit["NGIC"]
	
	# Load the config file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open config file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	scenario = my_config["rcp_scenario"]
	targyears = my_config["targyears"]
	baseyear = my_config["baseyear"]

	
	# Evenly sample an inverse normal distribution and permutate it
	# Note: This may be a bug being ported over from Kopp 2014 which could result in 
	# 		overconfident projections
	rng = np.random.default_rng(seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = np.full((T.shape[0], nsamps), np.nan)
	for i in np.arange(0,T.shape[0]):
		norm_inv_perm[i,:] = rng.permutation(norm_inv)
	
	## Generate the samples --------------------------------------------------------------
	# Initialize variable to hold the samples
	gicsamps = np.full((nsamps, T.shape[0], len(targyears)), np.nan)
	
	# Loop over the target years
	for i in np.arange(0,len(targyears)):
		
		# Values to use for producing samples at this target year
		thisYear = targyears[i]
		thisMeanGIC = meanGIC[:,i]
		thisT = T[:,:,i]
		thisNGIC = NGIC[i]
		
		# Generate the samples for this year
		if(thisNGIC > 0):
			temp = t.ppf(norm.cdf(norm_inv_perm), np.min((np.inf,thisNGIC-1))).T
			gicsamps[:,:,i] = np.dot(temp, thisT) + thisMeanGIC
		else:
			gicsamps[:,:,i] = np.nan
	
	# Reference these projections to the base year
	baseyear_idx = np.flatnonzero(targyears == baseyear)
	gicsamps = gicsamps - gicsamps[:,:,baseyear_idx]

	# Save the global glacier and ice caps projections to a pickle
	output = {"gicsamps": gicsamps}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Sum over all the regions
	total_glac_samps = np.apply_along_axis(np.sum, 1, gicsamps)
		
	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim  = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var  = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var  = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var  = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from glaciers and ice caps according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {}; scenario: {}; baseyear: {}".format(pipeline_id, scenario, baseyear)
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nsamps)
	samps[:,:,:] = total_glac_samps[:,:,np.newaxis]
	lat_var[:] 	 = np.inf
	lon_var[:] 	 = np.inf
	loc_var[:] 	 = -1

	# Close the netcdf
	rootgrp.close()	
	
	# Done
	return(0)



if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glacier projection stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', '-n', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', '-s', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	kopp14_project_glaciers(args.nsamps, args.seed, args.pipeline_id)
	
	exit()