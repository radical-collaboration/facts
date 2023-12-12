# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19
# Adapted for use in FACTS by Gregory Garner 20 November 2019

import os
import argparse
import pickle
from netCDF4 import Dataset
import numpy as np
from Import2lmData import *

class ProjectionError(Exception):
	pass
	
def endofhistory():
	return 2006

def ar5_preprocess_glaciers(scenario, startyr, tlm_flag, pipeline_id, climate_fname):
	
	# Define the input data directory
	indir = os.path.dirname(__file__)
	
	# Load the two-layer model data
	if tlm_flag:
		
		# Import the data
		tlm_dict = Import2lmData("surface_temperature", scenario, indir, refyear_start=1986, refyear_end=2005, climate_fname=climate_fname)
		
		# Filter the data for the appropriate years
		filtered_data_dict = Filter2lmData(tlm_dict, filter_years=np.arange(startyr,2301))
		
		# Extract the years
		data_years = filtered_data_dict['years']
		
		# FOR THE TEMPORARY 2LM DATA ONLY - USES ONLY 44 UNIQUE TRAJECTORIES
		# Find the unique temperature trajectories
		#temp_samples = np.unique(filtered_data_dict["samples"], axis=0)
		temp_samples = filtered_data_dict["samples"]
		
		# Find the mean and sd of the ensemble
		temp_mean = np.nanmean(temp_samples, axis=0)
		temp_sd = np.nanstd(temp_samples, axis=0)
		
	else:
	
		# Define the input data files
		temp_mean_filename = "{0}_temperature_mean.nc".format(scenario)
		temp_sd_filename = "{0}_temperature_sd.nc".format(scenario)
	
		# Extract the data ------------------------------------------
		nc = Dataset(os.path.join(indir,temp_mean_filename), 'r')
		temp_years = 2006 + nc.variables['bound_time'][:,1]/360
		year_idx = np.flatnonzero(temp_years >= startyr)
		temp_mean = nc.variables['air_temperature'][year_idx]
		nc.close()

		nc = Dataset(os.path.join(indir,temp_sd_filename), 'r')
		temp_sd = nc.variables['air_temperature'][year_idx]
		nc.close()
	
		# Set the output data years to the temperature years
		data_years = temp_years	

	# Integrate temperature to obtain K yr at ends of calendar years
	# Note - The original code I believe performs a cumulative sum of the standard
	# deviations, which is not correct.  Below I provide a fix to that bug as well as
	# a replication of the bug for diagnostic purposes.
	# Note - JG makes an assumption here so that standard deviations are able to be summed
	# over time.
	inttemp_mean = np.cumsum(temp_mean)
	#inttemp_sd = np.sqrt(np.cumsum(temp_sd**2))  # Fix the bug
	inttemp_sd = np.cumsum(temp_sd)  # Replicate the bug
	inttemp_samples = np.cumsum(temp_samples, axis=1)
		
	# Store preprocessed data in pickles
	output = {'temp_mean': temp_mean, 'temp_sd': temp_sd, 'inttemp_mean': inttemp_mean, \
				'inttemp_sd': inttemp_sd, 'data_years': data_years, 'startyr': startyr, \
				'scenario': scenario, 'temp_samples': temp_samples, 'inttemp_samples': inttemp_samples}
	
	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(0)


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glaciers pre-processing stage for the AR5 Global SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Scenario choices
	#scenario_choices = ['rcp85', 'rcp60', 'rcp45', 'rcp26', \
	#					'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp434', 'ssp460', 'ssp585']
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="Scenario [default=\'rcp85\']", default='rcp85')
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--baseyear', help="Year from which to start integrating temperature [default=2006]", type=int, default=2006)
	parser.add_argument('--tlm_data', help="Use the two-layer model data [default=1, use 2lm data]", default=1, type=int)
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ar5_preprocess_glaciers(args.scenario, args.baseyear, args.tlm_data, args.pipeline_id, args.climate_data_file)
	
	exit()