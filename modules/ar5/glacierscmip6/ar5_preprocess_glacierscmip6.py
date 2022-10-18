# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19
# Adapted for use in FACTS by Gregory Garner 20 November 2019

import os
import argparse
import pickle
import fnmatch
import re
from netCDF4 import Dataset
import numpy as np



class ProjectionError(Exception):
	pass
	
	
	
def endofhistory():
	return 2014
	
	
	
def FindInputModels(tasdir, scenario):
	
	# Acceptable SSP scenarios
	ssp_scenarios = ['ssp585', 'ssp370', 'ssp245', 'ssp126', 'ssp119']
	
	# Initialize models and scenarios to include
	include_models = []
	include_scenarios = []
	
	# Test the provided scenario
	scenario_test = re.search("^tlim(\d*\.?\d+)win(\d*\.?\d+)$", scenario)
	if(scenario_test):
		
		# This is a temperature limit, so extract the limit from the scenario string
		temp_target = float(scenario_test.group(1))
		temp_target_window = float(scenario_test.group(2))
		
		# Produce a list of models and scenarios that match the criteria
		(include_models, include_scenarios) = tas_limit_filter(tasdir, temp_target, temp_target_window)
		
	elif(scenario in ssp_scenarios):
		
		# Loop through the list of models available in the tas directory
		for this_model in os.listdir(tasdir):
			
			# Skip hidden directories and files
			if(re.search("^\.", this_model)):
				continue
				
			# Find the appropriate tas file for this model and scenario
			tas_filename = fnmatch.filter(os.listdir(os.path.join(tasdir, this_model)), "*{}*".format(scenario))[0]
			
			# If the file exists, append this model to the model list
			if(os.path.isfile(os.path.join(tasdir, this_model, tas_filename))):
				include_models.append(this_model)
				include_scenarios.append(scenario)
	
	else:
		
		# This is an invalid scenario
		raise Exception("Invalid scenario definition: {}".format(scenario))
	
	return(include_models, include_scenarios)



def tas_limit_filter(tasdir, temp_target, temp_target_window=0.25, ref_syear=1850, ref_eyear=1900, extrap_eyear=2099):
	
	# Initialize a running list of models and scenarios to include in the analysis
	include_models = []
	include_scenarios = []
	
	# Get a list of models available from the subdirectories available in the parent model directory
	modeldirs = os.listdir(tasdir)
	
	# Loop through the models
	for this_modeldir in modeldirs:
		
		# Skip any hidden directories or files
		if(re.search("^\.", this_modeldir)):
			continue
		
		# Open the historical file if available
		hist_filename = fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*historical*")[0]
		try:
			nc_hist = Dataset(os.path.join(tasdir, this_modeldir, hist_filename), 'r')
		except:
			continue
		
		# Get the years out of the historical file
		hist_years = nc_hist.variables['year'][:]
		
		# Determine which indices are in our reference average window
		ref_idx = np.flatnonzero(np.logical_and(hist_years >= ref_syear, hist_years <= ref_eyear))
		
		# Extract the temperature data from the reference period and produce the average
		ref_tas = np.mean(nc_hist.variables['tas'][ref_idx])
		
		# Close the historical netCDF file
		nc_hist.close()
		
		# Loop through all the remaining files in this model directory
		for this_filename in fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*ssp*"):
			
			# Open the netCDF file if available
			this_file = os.path.join(tasdir, this_modeldir, this_filename)
			nc = Dataset(this_file, 'r')
			
			# Get the years available in this projection
			proj_years = nc.variables['year'][:]
			
			# Determine the indices needed to calculate the rate of the 19-yr average
			stop_extrap_idx = np.flatnonzero(np.logical_and(proj_years >= (extrap_eyear-19), proj_years <= extrap_eyear))
			start_extrap_idx = np.flatnonzero(np.logical_and(proj_years >= (extrap_eyear-20), proj_years <= (extrap_eyear-1)))
			
			# Calculate the means and take the difference to get the rate
			start_extrap = np.mean(nc.variables['tas'][start_extrap_idx])
			stop_extrap = np.mean(nc.variables['tas'][stop_extrap_idx])
			tas_rate = stop_extrap - start_extrap
			
			# Extrapolate that to get a 19-yr average centered on 2100 and subtract the reference
			tas_mean_2100 = (stop_extrap + (2100 - extrap_eyear + 9) * tas_rate) - ref_tas
			
			# If tas_mean_2100 is within 0.25 deg of the temperature target, set this model as a keeper
			if(tas_mean_2100 >= temp_target - temp_target_window and tas_mean_2100 <= temp_target + temp_target_window):
				
				# Add the model to the list
				include_models.append(this_modeldir)
				
				# Add the scenario to the list
				this_scenario = re.search(r"(ssp\d{3})", this_filename).group(1)
				include_scenarios.append(this_scenario)
		
	return(include_models, include_scenarios)



def AggregateTAS(tasdir, modellist, scenariolist, ref_syear=1986, ref_eyear=2005, data_syear=2006, data_eyear=2100):
	
	# Make sure the reference years are within the historical data
	if(ref_eyear < ref_syear):
		raise Exception("ref_syear is later than ref_eyear")
	if(ref_syear < 1850 or ref_syear > 2014):
		raise Exception("ref_syear is beyond the limits of the historical data (1850 - 2014)")
	if(ref_eyear < 1850 or ref_eyear > 2014):
		raise Exception("ref_eyear is beyond the limits of the historical data (1850 - 2014)")
	
	# Make sure the projection years are within the projection data
	if(data_eyear < data_syear):
		raise Exception("data_syear is later than data_eyear")
	if(data_syear < 1850 or data_syear > 2099):
		raise Exception("data_syear is beyond the limits of the available data (1850 - 2099)")

	# Initialize the data structure to hold the temperature data
	model_temp = []
	
	# Loop through the models and scenarios
	for i in np.arange(len(modellist)):
		
		# Extract this model/scenario pair
		this_modeldir = modellist[i]
		this_scenario = scenariolist[i]
		
		# Open the historical file if available
		hist_filename = fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*historical*")[0]
		try:
			nc_hist = Dataset(os.path.join(tasdir, this_modeldir, hist_filename), 'r')
		except:
			raise Exception("Cannot open historical file for this model/scenario combination - {0} {1}".format(this_modeldir, this_scenario))
				
		# Get the years out of the historical file
		hist_years = nc_hist.variables['year'][:]
		
		# Determine which indices are in our reference average window
		ref_idx = np.flatnonzero(np.logical_and(hist_years >= ref_syear, hist_years <= ref_eyear))
		
		# Extract the temperature data from the reference period and produce the average
		ref_tas = np.mean(nc_hist.variables['tas'][ref_idx])
		
		# If the data start year is in the historical, extract the appropriate amount of data
		# to prepend to the projection data
		if(data_syear <= 2014):
			
			# Find the start of the data
			start_idx = np.flatnonzero(hist_years == data_syear)[0]
			
			# Extract the data
			hist_tas = nc_hist.variables['tas'][start_idx:]
			
		else:
			
			# There's no historical data to append for this data start year
			hist_tas = []
			
		
		# Close the historical netCDF file
		nc_hist.close()		
			
		# Open the netCDF file of the projection data
		proj_filename = fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*{}*".format(this_scenario))[0]
		proj_file = os.path.join(tasdir, this_modeldir, proj_filename)
		nc = Dataset(proj_file, 'r')
		
		# Get the years available in this projection
		proj_years = nc.variables['year'][:]
		
		# Define the index at which to start extracting the tas data
		if(data_syear >= 2015):
			start_idx = np.flatnonzero(proj_years == data_syear)[0]
		else:
			start_idx = 0
		
		# If the requested end year is greater than the maximum projected year, apply a
		# simple extrapolation to the requested end year
		# NOTE: This should be reasonable if extrapolating a couple of years, beyond
		#		that, the user should exercise caution!
		if(data_eyear > np.amax(proj_years)):
			
			# Find out how many years to extrapolate
			nextrapyears = data_eyear - np.amax(proj_years)
			
			# Read in the TAS data
			this_tas = nc.variables['tas'][start_idx:]
			
			# Append any historical data to this
			this_tas = np.append(hist_tas, this_tas)
			
			# Produce the annual means from the monthly data
			this_annual_tas = [np.mean(this_tas[x:x+12]) for x in np.arange(0, len(this_tas), 12)]
			
			# Linearly extrapolate to the end year
			trend = np.diff(this_annual_tas[-2:])
			extrap_trend = np.arange(nextrapyears+1)[1:] * trend
			extrap_tas = this_annual_tas[-1] + extrap_trend
			
			# Append the extrapolation to the annual tas values
			this_annual_tas = np.append(this_annual_tas, extrap_tas)
		
		else:
		
			# Determine the index that matches the requested end projection years
			stop_idx = np.flatnonzero(proj_years == data_eyear)[-1]
			
			# Read in the TAS data
			this_tas = nc.variables['tas'][start_idx:stop_idx]
			
			# Append any historical data to this
			this_tas = np.append(hist_tas, this_tas)
			
			# Produce the annual means from the monthly data
			this_annual_tas = [np.mean(this_tas[x:x+12]) for x in np.arange(0, len(this_tas), 12)]
		
		# Center the annual means
		this_annual_tas = this_annual_tas - ref_tas
		
		# Stack the annual temperatures from this model with the rest of the models
		model_temp.append(this_annual_tas)
		
		# Close this netCDF file
		nc.close()
	
	# Convert the model temperature data structure into a numpy array
	model_temp = np.array(model_temp)
	
	# Take the multi-model mean and standard deviation
	temp_mean = np.mean(model_temp, axis=0)
	temp_sd = np.std(model_temp, axis=0)
	
	# Provide an array of the years for which data is available
	datayears = np.arange(data_syear, data_eyear+1)
	
	# Done, return the quantities
	return(temp_mean, temp_sd, datayears)


def ar5_preprocess_glacierscmip6(scenario, startyr, pipeline_id):
	
	
	# Define the list of input models and scenarios
	tasdir = os.path.join(os.path.dirname(__file__), "tas")
	(include_models, include_scenarios) = FindInputModels(tasdir, scenario)
	
	if(not include_models):
		raise Exception("No models found for this scenario or temperature target - {}".format(scenario))
	
	# Find the mean and sd of the matched models/scenarios
	(temp_mean, temp_sd, data_years) = AggregateTAS(tasdir, include_models, include_scenarios, data_syear=startyr)

	# Integrate temperature to obtain K yr at ends of calendar years
	# Note - The original code I believe performs a cumulative sum of the standard
	# deviations, which is not correct.  Below I provide a fix to that bug as well as
	# a replication of the bug for diagnostic purposes.
	inttemp_mean = np.cumsum(temp_mean)
	#inttemp_sd = np.sqrt(np.cumsum(temp_sd**2))  # Assume independence across models
	inttemp_sd = np.cumsum(temp_sd)  # Assume correlation
		
	# Store preprocessed data in pickles
	output = {'temp_mean': temp_mean, 'temp_sd': temp_sd, 'inttemp_mean': inttemp_mean, \
				'inttemp_sd': inttemp_sd, 'data_years': data_years, 'startyr': startyr, \
				'scenario': scenario, 'include_models': include_models, 'include_scenarios': include_scenarios}
	
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
	parser.add_argument('--baseyear', help="Year from which to start integrating temperature [default=2006]", type=int, default=2006)
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ar5_preprocess_glacierscmip6(args.scenario, args.baseyear, args.pipeline_id)
	
	exit()