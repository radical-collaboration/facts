import numpy as np
import pickle
import sys
import os
import argparse
import time
from netCDF4 import Dataset
from scipy.stats import norm
from scipy.stats import t

''' kopp14_project_thermalexp.py

This script runs the thermal expansion projection task for the Kopp 2014 workflow. 
This task generates global contributions to sea-level change due to thermosteric heating
of the ocean.

Parameters: 
nsamps = Numer of samples to produce
seed = Seed for the random number generator
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickle file containing the global sea-level rise projections

'''

def kopp14_project_thermalexp(nsamps, seed, pipeline_id):
	
	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	rcp_scenario = my_config["rcp_scenario"]
	#datayears = my_config["datayears"]
	targyears = my_config["targyears"]
	#mergeZOSZOSTOGA = my_config["mergeZOSZOSTOGA"]
	#smoothwin = my_config["smoothwin"]
	#driftcorr = my_config["driftcorr"]
	#baseyear = my_config["baseyear"]
	GCMprobscale = my_config["GCMprobscale"]
	
	# Load the fit file
	fitfile = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(fitfile, 'rb')
	except:
		print("Cannot open fit file\n")
	
	# Extract the fit variables
	my_fit = pickle.load(f)
	f.close()
	
	ThermExpMean = my_fit["ThermExpMean"]
	ThermExpStd = my_fit["ThermExpStd"]
	ThermExpYears = my_fit["ThermExpYears"]	
	ThermExpN = my_fit["ThermExpN"]
	
	# Evenly sample an inverse normal distribution and permutate it
	# Note: This may be a bug being ported over from Kopp 2014 which could result in 
	# 		overconfident projections
	np.random.seed(seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = np.random.permutation(norm_inv)
	
	# Determine the scale coefficient
	ThermExpScale = norm.ppf(0.95)/norm.ppf(GCMprobscale)
	
	## Generate the samples --------------------------------------------------------------
	# Initialize variable to hold the samples
	thermsamps = np.full((nsamps, len(targyears)), np.nan)
	
	# Loop over the target years
	for i in np.arange(0,len(targyears)):
		
		# Find the index of ThermExp* that matches this target year
		this_year_ind = np.flatnonzero(ThermExpYears == targyears[i])
		
		# Generate the samples for this year
		#temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpN[this_year_ind]-1)
		temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpN[i]-1)  # Replicates bug in K14 master code
		thermsamps[:,i] = (ThermExpScale * temp * ThermExpStd[this_year_ind]) + ThermExpMean[this_year_ind]

	# Save the global thermal expansion projections to a pickle
	output = {"thermsamps": thermsamps, "targyears": targyears, "rcp_scenario": rcp_scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
	samp_dim = rootgrp.createDimension("samples", nsamps)

	# Populate dimension variables
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	samps = rootgrp.createVariable("samps", "f4", ("years", "samples"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from thermal expansion according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, rcp_scenario)
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nsamps)
	samps[:,:] = np.transpose(thermsamps)

	# Close the netcdf
	rootgrp.close()	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion projection stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', '-n', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', '-s', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	kopp14_project_thermalexp(args.nsamps, args.seed, args.pipeline_id)
	
	exit()