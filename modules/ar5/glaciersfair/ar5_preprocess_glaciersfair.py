# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19
# Adapted for use in FACTS by Gregory Garner 20 November 2019

import os
import argparse
import pickle
import fnmatch
import re
from netCDF4 import Dataset
from import_data import *
from filter_data import filter_data
from Smooth import Smooth
import numpy as np



class ProjectionError(Exception):
	pass
	
	
	
def endofhistory():
	return 2014
	


def tas_limit_filter(tasdict, temp_target, temp_target_window=0.25, ref_syear=1850, ref_eyear=1900):
	
	# Initialize a running list of models and scenarios to include in the analysis
	out_dict = {}
	
	# Determine which indices are in our reference average window
	ref_idx = np.flatnonzero(np.logical_and(tasdict["years"] >= ref_syear, tasdict["years"] <= ref_eyear))
	
	# Extract the temperature data from the reference period and produce the average
	ref_tas = np.mean(tasdict['data'][:,ref_idx], axis=1)
	
	# Smooth the temperature data over the 19 year window
	tas_smoothed = np.apply_along_axis(Smooth, axis=1, arr=tasdict["data"], w=19)
	
	# Find the tas value to test for the temperature target
	tas_test = tas_smoothed[:,-1]
	#tas_test = np.nanmax(tas_smoothed, axis=1)

	# Get the indices where the 2100 temperature falls within the temperature limit window
	match_idx = np.flatnonzero(np.logical_and(tas_test >= temp_target - temp_target_window, tas_test <= temp_target + temp_target_window))
	
	# Subset the original tasdict for the matched indices
	out_dict["ensemble"] = tasdict["ensemble"][match_idx]
	out_dict["GCM"] = tasdict["GCM"][match_idx]
	out_dict["scenario"] = tasdict["scenario"][match_idx]
	out_dict["years"] = tasdict["years"]
	out_dict["data"] = tasdict["data"][match_idx,:]

	# Return the output dictionary
	return(out_dict)




def ar5_preprocess_glaciersfair(scenario, startyr, pipeline_id):
	
	# Define the temperature input file name
	infilename = "CLIMATE_FORCING_1850.csv"
	infile = os.path.join(os.path.dirname(__file__), infilename)
	
	# Acceptable SSP scenarios
	ssp_scenarios = ['ssp585', 'ssp370', 'ssp245', 'ssp126', 'ssp119']
	
	# Import the temperature data
	temp_data = import_data(infile)
	
	# Test the provided scenario
	scenario_test = re.search("^tlim(\d*\.?\d+)win(\d*\.?\d+)$", scenario)
	if(scenario_test):
		
		# This is a temperature limit, so extract the limit from the scenario string
		temp_target = float(scenario_test.group(1))
		temp_target_window = float(scenario_test.group(2))
		
		# Produce a list of models and scenarios that match the criteria
		temp_data_filtered = tas_limit_filter(temp_data, temp_target, temp_target_window)
		
	elif(scenario in ssp_scenarios):
		
		# Filter the temperature data for this particular scenario
		temp_data_filtered = filter_data(temp_data, ensemble="FAIR", scenario=scenario.upper())
	
	else:
		
		# This is an invalid scenario
		raise Exception("Invalid scenario definition: {}".format(scenario))
	
	# The module is calibrated to use the temperature reference period for AR5, so center
	# the temperature data to the mean of that period
	ref_idx = np.flatnonzero(np.logical_and(temp_data_filtered["years"] >= 1986, temp_data_filtered["years"] <= 2005))
	ref_tas = np.nanmean(temp_data_filtered["data"][:,ref_idx], axis=1)
	temp_data_filtered["data"] = temp_data_filtered["data"] - ref_tas[:,np.newaxis]
	
	# Find the mean and sd of the matched models/scenarios
	temp_mean = np.nanmean(temp_data_filtered['data'], axis=0)
	temp_sd = np.nanstd(temp_data_filtered['data'], axis=0)
	data_years = temp_data_filtered['years']
	
	# Find which year in the data years is the start year
	baseyear_idx = np.flatnonzero(data_years == startyr)

	# Integrate temperature to obtain K yr at ends of calendar years
	# Note - The original code I believe performs a cumulative sum of the standard
	# deviations, which is not correct.  Below I provide a fix to that bug as well as
	# a replication of the bug for diagnostic purposes.
	inttemp_mean = np.cumsum(temp_mean)
	#inttemp_sd = np.sqrt(np.cumsum(temp_sd**2))  # Assume independence across models
	inttemp_sd = np.cumsum(temp_sd)  # Assume correlation
	
	# Integrated quantities must be centered on the baseline year
	inttemp_mean -= inttemp_mean[baseyear_idx]
	inttemp_sd -= inttemp_sd[baseyear_idx]
		
	# Store preprocessed data in pickles
	output = {'temp_mean': temp_mean, 'temp_sd': temp_sd, 'inttemp_mean': inttemp_mean, \
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
	parser = argparse.ArgumentParser(description="Run the glaciers pre-processing stage for the AR5 Global SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="SSP scenario (i.e. ssp585) or temperature target (i.e. tlim2.0win0.25)", default='ssp585')
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--startyear', help="Year from which to start integrating temperature [default=2005]", type=int, default=2005)
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ar5_preprocess_glaciersfair(args.scenario, args.startyear, args.pipeline_id)
	
	exit()