import numpy as np
import argparse
import pickle
import os
import sys
import copy
from netCDF4 import Dataset

''' dp_preprocess_icesheet.py

Preprocess data for use in the dp16 icesheet module.

Parameters:
scenario - Emissions scenario of interest

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def dp16_preprocess_icesheet(scenario, baseyear, pipeline_id):

	# Define the input file names based on scenario
	eais_filename = "dp16_eais_{}.nc".format(scenario)
	wais_filename = "dp16_wais_{}.nc".format(scenario)

	# Get the years
	years = LoadNetCDF(eais_filename, "years")

	# Get the actual data
	eais_samps = LoadNetCDF(eais_filename, "samps")
	wais_samps = LoadNetCDF(wais_filename, "samps")

	# Get the values for the baseyear of interest
	eais_refs = np.apply_along_axis(FindRefVals, axis=0, arr=eais_samps, years=years, baseyear=baseyear)
	wais_refs = np.apply_along_axis(FindRefVals, axis=0, arr=wais_samps, years=years, baseyear=baseyear)

	# Center the samples to the base year
	eais_samps -= eais_refs
	wais_samps -= wais_refs

    ###################################################
    # Store the data in a pickle
	output = {'years': years, 'eais_samps': eais_samps, 'wais_samps': wais_samps, \
				'scenario': scenario, 'baseyear': baseyear}

	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfilename = "{}_data.pkl".format(pipeline_id)
	outfile = open(os.path.join(outdir, outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()


def FindRefVals(timeseries, years, baseyear):

	# Append a zero to the beginning of the timeseries at year 2000
	timeseries = np.append(np.array([0.0]), timeseries)
	years = np.append(np.array([2000]), years)

	# Interpolate to the appropriate base year
	ref_val = np.interp(baseyear, years, timeseries, left=0.0)

	# Return the value
	return(ref_val)



def LoadNetCDF(filename, variable):

	# Open the file
	nc = Dataset(filename, 'r')

	# Extract the variable
	var = nc.variables[variable][...]

	# Close the file
	nc.close()

	# Return the variable
	return(var)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the DP16 icesheet pre-processing stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="Emissions scenario of interest [default = rcp85]", choices=["rcp85", "rcp45", "rcp26"], default="rcp85")
	parser.add_argument('--baseyear', help="Year to which projections are referenced [default = 2000]", default=2000, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# If baseyear is before year 2000, fail because we have no data prior to year 2000
	if args.baseyear < 2000:
		raise Exception("Base year must be greater than 2000: baseyear = {}".format(args.baseyear))

	# Run the preprocessing stage with the provided arguments
	dp16_preprocess_icesheet(args.scenario, args.baseyear, args.pipeline_id)

	exit()
