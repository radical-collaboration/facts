import numpy as np
import os
import sys
import re
import fnmatch
from netCDF4 import Dataset

'''
tas_limit_filter()

Returns a list of models and a list of corresponding scenarios to be used to filter the
ocean dynamics module for a specific temperature target.
'''

def tas_limit_filter(tasdir, temp_target, temp_target_window=0.25, ref_syear=1875, ref_eyear=1900, extrap_eyear=2099):
	
	# Initialize a running list of models and scenarios to include in the analysis
	include_models = []
	include_scenarios = []
	temp_anom = []
	
	# Get a list of models available from the subdirectories available in the parent model directory
	modeldirs = os.listdir(tasdir)
	
	# Loop through the models
	for this_modeldir in modeldirs:
		
		# Next if this is a hidden file or directory
		if(re.search(r"^\.", this_modeldir)):
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
			tas_mean_2100 = (stop_extrap + (2100 - extrap_eyear - 9) * tas_rate) - ref_tas
			
			# If tas_mean_2100 is within 0.25 deg of the temperature target, set this model as a keeper
			if(tas_mean_2100 >= temp_target - temp_target_window and tas_mean_2100 <= temp_target + temp_target_window):
				
				# Add the model to the list
				include_models.append(this_modeldir)
				
				# Add the scenario to the list
				this_scenario = re.search(r"(ssp\d{3})", this_filename).group(1)
				include_scenarios.append(this_scenario)
				
				# Add the temperature to the list
				temp_anom.append(tas_mean_2100)
		
	return(include_models, include_scenarios, temp_anom)


if __name__ == "__main__":
	
	dir = "/Users/ggg46/research/slr_framework/code/facts/modules/ipccar6/oceandynamics/cmip6/tas"
	
	temp_target = 5.0
	temp_target_window = 0.5
	ref_syear = 1850
	ref_eyear = 1900
	
	(models, scenarios, temp) = tas_limit_filter(dir, temp_target, temp_target_window, ref_syear=ref_syear, ref_eyear=ref_eyear)
	
	print("Limit = {0}  Window = {1}  Baseline Years = {2}-{3}".format(temp_target, temp_target_window, ref_syear, ref_eyear))
	
	for i in np.arange(len(models)):
		
		print("  | {0} - {1} : {2:.3f}".format(models[i], scenarios[i], temp[i]))
	
	print("Mean temperature anomaly: {0:.3f}".format(np.mean(temp)))
	
	exit()