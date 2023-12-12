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


''' TLM thermal expansion projection '''
def tlm_project_thermalexpansion(seed, nsamps, pipeline_id,scenario,pyear_start,pyear_end,pyear_step,baseyear):

	'''
	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")

	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	'''

	#scenario = my_config["scenario"]
	# targyears = my_config["targyears"]
	# baseyear = my_config["baseyear"]
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)


	# Set the RNG seed
	rng = np.random.default_rng(seed)

	# Load the preprocessed data
	data_file = "{}_tlmdata.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")

	# Extract the data variables
	my_data = pickle.load(f)
	f.close()

	ohc_samps = my_data['ohc_samps']
	scenario = my_data['scenario']
	data_years = my_data['data_years']

	# Load the Fitted data
	data_file = "{}_tlmfit.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")

	# Extract the data variables
	my_data = pickle.load(f)
	f.close()

	mean_expcoefs = my_data['mean_expcoefs']
	std_expcoefs = my_data['std_expcoefs']
	include_models = my_data['include_models']

	# Generate indices that sample from OHC values
	#ohc_samps_idx = rng.choice(np.arange(ohc_samps.shape[0]), nsamps)
	ohc_samps_idx = np.arange(ohc_samps.shape[0])
	ohc_samps = ohc_samps[ohc_samps_idx,:]

	# Generate samples assuming normal distribution
	expcoef_samps = rng.normal(loc=mean_expcoefs, scale=std_expcoefs, size=(nsamps,1))

	# Produce the projection samples
	gte_samps = ohc_samps * expcoef_samps

	# Center these samples on the baseyear
	baseyear_idx = np.flatnonzero(data_years == baseyear)
	gte_samps = gte_samps - gte_samps[:,baseyear_idx]

	# Subset the samples for the projection years
	targyear_idx = np.isin(data_years, targyears)
	gte_samps = gte_samps[:,targyear_idx]

	# Invert the dimensions of the variable and convert from m to mm
	gte_samps *= 1000.

	# Save the projections to a pickle
	output = {"thermsamps": gte_samps, "targyears": targyears, "baseyear": baseyear, \
				"include_models": include_models, "scenario": scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	nyr = len(targyears)
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
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = "Global SLR contribution from Thermal Expansion according to Two-Layer Model workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}".format(pipeline_id)
	rootgrp.scenario = scenario
	rootgrp.baseyear = baseyear
	rootgrp.comment = "Included Models: " + ",".join([str(x) for x in include_models])
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = gte_samps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the IPCC AR6 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)",
						default='ssp585')
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2000)

	# Parse the arguments
	args = parser.parse_args()

	# Intermediate file produced in this stage
	projfile = os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(args.pipeline_id))

	# If the intermediate file exists for this pipeline id, skip the projection stage
	if os.path.isfile(projfile):
		print("{} found, skipping projection stage".format(projfile))
		sys.exit()

	# Run the project stage with the user defined scenario
	tlm_project_thermalexpansion(args.seed, args.nsamps, args.pipeline_id, args.scenario, args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear)

	# Done
	sys.exit()