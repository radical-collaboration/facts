import os
import sys
import re
import numpy as np
from import_data import import_data
from filter_data import filter_data
from import_temp_data import import_temp_data
from filter_temp_data import filter_temp_data
from Smooth import Smooth


def FindFAIRInputSamples(forcing_dict, scenario):

	# Acceptable SSP scenarios
	ssp_scenarios = ['ssp585', 'ssp370', 'ssp245', 'ssp126', 'ssp119']

	# Test the provided scenario
	scenario_test = re.search("^tlim(\d*\.?\d+)win(\d*\.?\d+)$", scenario)
	if(scenario_test):

		# This is a temperature limit, so extract the limit from the scenario string
		temp_target = float(scenario_test.group(1))
		temp_target_window = float(scenario_test.group(2))

		# Produce a list of models and scenarios that match the criteria
		sample_dict = tas_limit_filter(forcing_dict, temp_target, temp_target_window)

	elif(scenario in ssp_scenarios):

		# Populate sample_dict with this scenario and the sample numbers
		n_unique_samples = len(np.unique(forcing_dict["GCM"]))
		sample_dict = {scenario: np.arange(1,n_unique_samples+1)}

	else:

		# This is an invalid scenario
		raise Exception("Invalid scenario definition: {}".format(scenario))

	return(sample_dict)



def tas_limit_filter(tasdict, temp_target, temp_target_window=0.25, ref_syear=1850, ref_eyear=1900):

	# Initialize a running list of models and scenarios to include in the analysis
	sample_dict = {}

	# Determine which indices are in our reference average window
	ref_idx = np.flatnonzero(np.logical_and(tasdict["years"] >= ref_syear, tasdict["years"] <= ref_eyear))

	# If we can't build a reference period, warn the user and use the raw data
	if len(ref_idx > 0):

		# Extract the temperature data from the reference period and produce the average
		ref_tas = np.mean(tasdict['data'][:,ref_idx], axis=1)

		# Subtract the reference period from the data
		tasdict['data'] = tasdict['data'] - ref_tas[:,None]

	else:

		print("Unable to find reference period in forcing data.  Using raw data.")

	# Smooth the temperature data over the 19 year window
	tas_smoothed = np.apply_along_axis(Smooth, axis=1, arr=tasdict["data"], w=19)

	# Find the tas value to test for the temperature target
	#tas_test = tas_smoothed[:,-1]
	tas_test = np.nanmax(tas_smoothed, axis=1)

	# Get the indices where the 2100 temperature falls within the temperature limit window
	match_idx = np.flatnonzero(np.logical_and(tas_test >= temp_target - temp_target_window, tas_test <= temp_target + temp_target_window))

	# Separate the matches by SSP scenario
	matched_ssps = np.unique(tasdict["scenario"][match_idx])
	for this_scenario in matched_ssps:

		# Indices for this scenario (sub-indices of match_idx)
		ssp_idx = np.flatnonzero(tasdict["scenario"][match_idx] == this_scenario)

		# Extract the samples for this particular scenario
		# Note: We're grabbing the FAIR ensemble member number from the entry in the GCM column (FAIR_XXX),
		# so the index 5 starts at the digit in the entry
		matched_sample = np.array([int(x[5:]) for x in tasdict["GCM"][match_idx[ssp_idx]]])

		# Store these samples in the output dictionary
		sample_dict[this_scenario.lower()] = matched_sample

	return(sample_dict)



def SubsetSLEProjections(sample_dict):

	# SLE projection directory
	#sle_dir = "191220_emulated"
	#sle_dir = "201011_proj_TIMESERIES"
	sle_dir = "2lm_projections"

	# Initialize the output dictionary
	sle_dict = {"ice_source": [], "region": [], "year": [], "scenario-sample": [],\
				"GSAT": [], "SLE": []}

	# Loop over the required SSPs from the sample dictionary
	for this_scenario in sample_dict.keys():

		# Open this matched scenario file
		filename = os.path.join(sle_dir, "projections_FAIR_{0}.csv".format(this_scenario.upper()))
		this_sle_dict = import_data(filename, "FAIR")

		# Filter this data for the appropriate samples
		this_sle_dict = filter_data(this_sle_dict, "FAIR", sample=sample_dict[this_scenario], ice_source="Glaciers")

		# Append these data to the output structure
		sle_dict["ice_source"].extend(this_sle_dict["ice_source"])
		sle_dict["region"].extend(this_sle_dict["region"])
		sle_dict["year"].extend(this_sle_dict["year"])
		sle_dict["GSAT"].extend(this_sle_dict["GSAT"])
		sle_dict["SLE"].extend(this_sle_dict["SLE"])

		# Add a field called "scenario-sample" to the dictionary
		scenario_sample = ["{0}-{1}".format(this_scenario, x) for x in this_sle_dict["sample"]]
		sle_dict["scenario-sample"].extend(scenario_sample)

	# Convert everything over into numpy arrays
	sle_dict["ice_source"] = np.array(sle_dict["ice_source"])
	sle_dict["region"] = np.array(sle_dict["region"])
	sle_dict["year"] = np.array(sle_dict["year"])
	sle_dict["GSAT"] = np.array(sle_dict["GSAT"])
	sle_dict["SLE"] = np.array(sle_dict["SLE"])
	sle_dict["scenario-sample"] = np.array(sle_dict["scenario-sample"])

	# Return the sea level projection dictionary
	return(sle_dict)


if __name__ == "__main__":

	#filename = "CLIMATE_FORCING_1850.csv"
	filename = "20201009_CLIMATE_FORCING.csv"

	x = import_temp_data(filename)

	x_filtered = filter_temp_data(x, ensemble="FAIR", scenario=["SSP119", "SSP126", "SSP245", "SSP370", "SSP585"])

	sample_dict = FindFAIRInputSamples(x_filtered, "tlim2.0win0.25")
	#sample_dict = FindFAIRInputSamples(x_filtered, "ssp585")

	print(sample_dict["ssp126"])


	#sle_dict = SubsetSLEProjections(sample_dict)

	exit()
