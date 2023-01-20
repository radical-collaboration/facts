# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19
# Adapted for use in FACTS by Gregory Garner 20 November 2019

import os
import argparse
import pickle
from netCDF4 import Dataset
import numpy as np

class ProjectionError(Exception):
	pass
	
def endofhistory():
	return 2006

def ar5_preprocess_thermalexpansion(scenario, startyr, pipeline_id):
	
	# Define the input data files
	indir = os.path.dirname(__file__)
	temp_mean_filename = "{0}_temperature_mean.nc".format(scenario)
	temp_sd_filename = "{0}_temperature_sd.nc".format(scenario)
	exp_mean_filename = "{0}_expansion_mean.nc".format(scenario)
	exp_sd_filename = "{0}_expansion_sd.nc".format(scenario)
	
	# Extract the data ------------------------------------------
	nc = Dataset(os.path.join(indir,exp_mean_filename), 'r')
	exp_years = 2007 + nc.variables['time'][:]/360  # Days since 2007-1-1
	year_idx = np.flatnonzero(exp_years >= startyr)	
	exp_mean = nc.variables['global_average_thermosteric_sea_level_change'][year_idx]
	nc.close()
	
	nc = Dataset(os.path.join(indir,exp_sd_filename), 'r')
	exp_sd = nc.variables['global_average_thermosteric_sea_level_change'][year_idx]
	nc.close()
	
	nc = Dataset(os.path.join(indir,temp_mean_filename), 'r')
	temp_years = 2006 + nc.variables['bound_time'][:,1]/360
	year_idx = np.flatnonzero(temp_years >= startyr)
	temp_mean = nc.variables['air_temperature'][year_idx]
	nc.close()

	nc = Dataset(os.path.join(indir,temp_sd_filename), 'r')
	temp_sd = nc.variables['air_temperature'][year_idx]
	nc.close()
	
	# Set the output data years to the expansion years
	data_years = exp_years	

	# Integrate temperature to obtain K yr at ends of calendar years
	# Note - The original code I believe performs a cumulative sum of the standard
	# deviations, which is not correct.  Below I provide a fix to that bug as well as
	# a replication of the bug for diagnostic purposes.
	inttemp_mean = np.cumsum(temp_mean)
	#inttemp_sd = np.sqrt(np.cumsum(temp_sd**2))  # Fix the bug
	inttemp_sd = np.cumsum(temp_sd)  # Replicate the bug
		
	# Store preprocessed data in pickles
	output = {'temp_mean': temp_mean, 'temp_sd': temp_sd, 'exp_mean': exp_mean,\
				'exp_sd': exp_sd, 'inttemp_mean': inttemp_mean, \
				'inttemp_sd': inttemp_sd, 'data_years': data_years, 'startyr': startyr, \
				'scenario': scenario}
	
	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(0)


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion pre-processing stage for the AR5 Global SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Scenario choices
	#scenario_choices = ['rcp85', 'rcp60', 'rcp45', 'rcp26', 'sresa1b', \
	#					'ssp126', 'ssp129', 'ssp245', 'ssp370', 'ssp434', 'ssp460', 'ssp585']
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="Scenario [default=\'rcp85\']",  default='rcp85')
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--baseyear', help="Year from which to start integrating temperature [default=2006]", type=int, default=2006)
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ar5_preprocess_thermalexpansion(args.scenario, args.baseyear, args.pipeline_id)
	
	exit()