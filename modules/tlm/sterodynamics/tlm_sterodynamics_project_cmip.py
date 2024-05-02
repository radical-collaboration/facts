import numpy as np
import pickle
import sys
import os
import argparse
import time
from netCDF4 import Dataset
from scipy.stats import norm
from scipy.stats import t

''' tlm_project_oceandynamics.py

This runs the projection stage for the ocean dynamics component of the IPCC AR6
workflow. The global projections are composed only of the thermal expansion portion.
Localized projections for ocean dynamics are done in post-processing to remain consistent
with the other module in this module set.

Parameters:
nsamps = Numer of samples to produce
seed = Seed for the random number generator
pipeline_id = Unique identifier for the pipeline running this code

Note that the value of 'nsamps' and 'seed' are passed to both the projection stage and
post-processing stage when run within FACTS.
'''

def tlm_project_cmip(nsamps, seed, pipeline_id):

	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")

	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()

	scenario = my_config["scenario"]
	targyears = my_config["targyears"]
	baseyear = my_config["baseyear"]
	GCMprobscale = my_config["GCMprobscale"]

	# Load the ZOSTOGA file
	configfile = "{}_ZOSTOGA.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open ZOSTOGA file\n")

	# Extract the ZOSTOGA variables
	my_data = pickle.load(f)
	f.close()

	zostoga_modellist = my_data["zostoga_modellist"]
	zostoga_scenariolist = my_data["zostoga_scenariolist"]

	# Load the fit file
	fitfile = "{}_thermalexp_fit.pkl".format(pipeline_id)
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
	ThermExpDOF = my_fit["ThermExpDOF"]

	# Evenly sample an inverse normal distribution and permutate it
	# Note: This may be a bug being ported over from Kopp 2014 which could result in
	# 		overconfident projections
	rng = np.random.default_rng(seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = rng.permutation(norm_inv)

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
		temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpDOF[this_year_ind])
		#temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpDOF[i])  # Replicates bug in K14 master code
		thermsamps[:,i] = (ThermExpScale * temp * ThermExpStd[this_year_ind]) + ThermExpMean[this_year_ind]

	# Save the global thermal expansion projections to a pickle
	output = {"thermsamps": thermsamps, "targyears": targyears, "scenario": scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=4)
	outfile.close()

	# Produce the included model string
	model_string_pieces = ["{0}-{1}".format(zostoga_modellist[x], zostoga_scenariolist[x]) for x in np.arange(len(zostoga_modellist))]
	model_string = "Models and scenarios included: " + ", ".join(model_string_pieces)

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	teyear_dim = rootgrp.createDimension("ThermExpYears", len(ThermExpYears))
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	teyear_var = rootgrp.createVariable("ThermExpYears", "i4", ("ThermExpYears",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)
	teexpmean = rootgrp.createVariable("ThermExpMean", "f4", ("ThermExpYears",))
	teexpstd = rootgrp.createVariable("ThermExpStd", "f4", ("ThermExpYears",))
	teexpn = rootgrp.createVariable("ThermExpN", "i4", ("ThermExpYears",))
	teexpdof = rootgrp.createVariable("ThermExpDOF", "i4", ("ThermExpYears",))

	# Assign attributes
	rootgrp.description = "Global SLR contribution from thermal expansion according to IPCC AR6 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}. ".format(pipeline_id) + model_string
	rootgrp.baseyear = baseyear
	rootgrp.scenario = scenario
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	teyear_var[:] = ThermExpYears
	samp_var[:] = np.arange(nsamps)
	teexpmean[:] = ThermExpMean
	teexpstd[:] = ThermExpStd
	teexpn[:] = ThermExpN
	teexpdof[:] = ThermExpDOF
	samps[:,:,:] = thermsamps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the IPCC AR6 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Intermediate file produced in this stage
	projfile = os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(args.pipeline_id))

	# If the intermediate file exists for this pipeline id, skip the projection stage
	if os.path.isfile(projfile):
		print("{} found, skipping projection stage".format(projfile))
		sys.exit()

	# Run the project stage with the user defined scenario
	tlm_project_cmip(args.nsamps, args.seed, args.pipeline_id)

	# Done
	sys.exit()