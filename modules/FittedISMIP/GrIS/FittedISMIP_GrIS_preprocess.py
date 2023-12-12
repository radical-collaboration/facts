import numpy as np
import argparse
import pickle
import os
import sys
from import_temp_data import *
from filter_temp_data import *
from Import2lmData import *

''' FittedISMIP_preprocess_icesheet.py

Preprocess data for use in the FittedISMIP icesheet module.

Parameters:
scenario - Emissions scenario of interest

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def FittedISMIP_preprocess_icesheet(scenario, tlm_flag, pipeline_id, climate_fname):

	# Load the two-layer model data?
	if tlm_flag != 0:
		# Load the data
		tlm_dict = Import2lmData(variable="surface_temperature", scenario=scenario, climate_fname=climate_fname)

		# Filter the data if necessary

		# Extract the pertinent fields to pass on to other stages
		years = tlm_dict['years']
		temp_data = tlm_dict['samples']

	else:
		# Temperature file name
		tempfile = os.path.join(os.path.dirname(__file__), "20201009_CLIMATE_FORCING.csv")

		# Load and filter the temperature data for this scenario
		temp_dict = import_temp_data(tempfile)
		filtered_temp_data = filter_temp_data(temp_dict, years=temp_dict['years'], scenario=scenario.upper())

		# Extract the pertinent fields
		years = filtered_temp_data['years']
		temp_data = filtered_temp_data['data']

    ###################################################
    # Store the data in a pickle
	output = {'years': years, 'temp_data': temp_data, 'scenario': scenario}

	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfilename = "{}_data.pkl".format(pipeline_id)
	outfile = open(os.path.join(outdir, outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the FittedISMIP icesheet pre-processing stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Scenario choices
	scenario_choices = ["rcp85", "rcp60", "rcp45", "rcp26", \
						"ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp585"]

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="Emissions scenario of interest [default = ssp585]", default="ssp585")
	parser.add_argument('--tlm_data', help="Use two-layer model temperature trajectories [default = 1, do not use]", choices=[0,1], default=1, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the provided arguments
	FittedISMIP_preprocess_icesheet(args.scenario, args.tlm_data, args.pipeline_id, args.climate_data_file)

	exit()
