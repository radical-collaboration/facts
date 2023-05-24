import argparse
import os
import sys
import re
import pickle as p
import numpy as np
from import_data import import_data
from filter_data import filter_data
from FindFAIRInputSamples import *

''' ipccar6_preprocess_ismipemuicesheet.py

This runs the preprocessing stage for the ISMIP Emulated ice sheet component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code
scenario = RCP scenario of interest

'''

def ipccar6_preprocess_ismipemuicesheet(pipeline_id, scenario, pyear_start, pyear_end, pyear_step, model_driver, baseyear):

	# Get the searchable RCP scenario format
	#scenario_dict = {'ssp585': "SSP585", 'ssp370': "SSP370", 'ssp245': "SSP245", 'ssp126': "SSP126", 'ssp119': "SSP119"}
	#scenario_search = scenario_dict[scenario]

	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Load the data
	if(model_driver == "CMIP6"):
		filename = os.path.join(os.path.dirname(__file__), "191208_annual_CMIP6", "projections_{}.csv".format(scenario.upper()))
		data_dict = import_data(filename, model_driver)
	else:
		filename = os.path.join(os.path.dirname(__file__), "20210215_CLIMATE_FORCING_IPCC.csv")
		forcing_data = import_temp_data(filename)
		forcing_data_filtered = filter_temp_data(forcing_data, ensemble="FAIR", scenario=["SSP119", "SSP126", "SSP245", "SSP370", "SSP585"])
		sample_dict = FindFAIRInputSamples(forcing_data_filtered, scenario)
		data_dict = SubsetSLEProjections(sample_dict)

	# Filter the input data for ice sheet and target years
	eais_dict = filter_data(data_dict, model_driver, ice_source="AIS", region="EAIS", year=targyears)
	wais_dict = filter_data(data_dict, model_driver, ice_source="AIS", region="WAIS", year=targyears)
	pen_dict = filter_data(data_dict, model_driver, ice_source="AIS", region="PEN", year=targyears)
	gis_dict = filter_data(data_dict, model_driver, ice_source="GrIS", region="ALL", year=targyears)

	# Generate the sample data structures
	eais_samps = MakeDataStructure(eais_dict, model_driver)
	wais_samps = MakeDataStructure(wais_dict, model_driver)
	pen_samps = MakeDataStructure(pen_dict, model_driver)
	gis_samps = MakeDataStructure(gis_dict, model_driver)

	# Define the linear trend terms to be added to the samples
	#trend_mean = {"EAIS": -0.02, "WAIS": 0.28, "PEN": 0.06, "GIS": 0.46}	# SOD
	#trend_sd = {"EAIS": 0.05, "WAIS": 0.03, "PEN": 0.01, "GIS": 0.04}		# SOD
	trend_mean = {"EAIS": 0.09, "WAIS": 0.18, "PEN": 0.06, "GIS": 0.19}
	trend_sd = {"EAIS": 0.04, "WAIS": 0.09, "PEN": 0.03, "GIS": 0.1}

	# Populate the output dictionary
	outdata = {'eais_samps': eais_samps, 'wais_samps': wais_samps, 'pen_samps': pen_samps, \
				'gis_samps': gis_samps, 'scenario': scenario, 'targyears': targyears, \
				'trend_mean': trend_mean, 'trend_sd': trend_sd, 'baseyear': baseyear, \
				'model_driver': model_driver}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()


def MakeDataStructure(data_dict, model_driver):

	# Initialize the data structure
	outdata = []

	# If this is CMIP6 data
	if(model_driver == "CMIP6"):

		# Get a list of the GCM, param_n, and emulator_n available
		all_gcm = np.unique(data_dict['GCM'])
		all_param_n = np.unique(data_dict['param_n'])
		all_emulator_n = np.unique(data_dict['emulator_n'])

		# Loop over the available GCMs, param_ns, and emualtor_ns
		for this_param_n in all_param_n:
			for this_emulator_n in all_emulator_n:

				# Get the available years
				filtered_dict = filter_data(data_dict, model_driver, param_n=this_param_n, emulator_n=this_emulator_n)
				all_years = np.unique(filtered_dict['year'])

				# Get the order of the years
				year_order = np.argsort(all_years)

				# Get the time series for this model
				ts_data = filtered_dict['SLE'][year_order]

				# Append this time series to the output data structure
				outdata.append(ts_data)

	# Otherwise, this is FAIR data
	else:

		# Get a list of the unique samples
		all_samples = np.unique(data_dict['scenario-sample'])

		# Loop over the available samples
		for this_sample in all_samples:

			# Get the available years
			filtered_dict = filter_data(data_dict, model_driver, scenario_sample = this_sample)
			all_years = np.unique(filtered_dict['year'])

			# Get the order of the years
			year_order = np.argsort(all_years)

			# Get the time series for this model
			ts_data = filtered_dict['SLE'][year_order]

			# Append this time series to the output data structure
			outdata.append(ts_data)

	# Make this a numpy array
	outdata = np.array(outdata)

	# Return the data structure
	return(outdata)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the IPCC AR6 ISMIP Emulated ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the IPCC AR6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--model_driver', help="The temperature model that drives the projections", default='FAIR', choices=['CMIP6', 'FAIR'])
	parser.add_argument('--scenario', help="SSP scenario of interest", default='ssp585')
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Baseline year to which contributions are zeroed", default=2005, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	ipccar6_preprocess_ismipemuicesheet(args.pipeline_id, args.scenario, args.pyear_start, args.pyear_end, args.pyear_step, args.model_driver, args.baseyear)

	# Done
	exit()
