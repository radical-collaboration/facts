import numpy as np
import argparse
import pickle
import os
import sys
import copy
from netCDF4 import Dataset

''' dp_preprocess_icesheet.py

Preprocess data for use in the dp21 icesheet module.

Parameters:
scenario - Emissions scenario of interest

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def dp21_preprocess_icesheet(scenario, baseyear, pipeline_id, climate_data_file = ''):

	if len(climate_data_file) > 0:
		scens=["rcp26","rcp45","rcp85"]
		years, eais_samps, wais_samps = ReadScenarioFile(scens[0],baseyear)
		eais_samps=eais_samps[:,:,np.newaxis]
		wais_samps=wais_samps[:,:,np.newaxis]

		for ii in range(1,len(scens)):
			years, e, w = ReadScenarioFile(scens[ii],baseyear)
			eais_samps=np.append(eais_samps,e[:,:,np.newaxis],axis=2)
			wais_samps=np.append(wais_samps,w[:,:,np.newaxis],axis=2)
	else:
		years, eais_samps, wais_samps = ReadScenarioFile(scenario,baseyear)

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

def ReadScenarioFile(scenario,baseyear):

		# Dictionary for mapping scenario names
		scen_dict = {"rcp85": "rcp85",
					"rcp45": "rcp45",
					"rcp26": "rcp26",
					"ssp585": "rcp85",
					"ssp245": "rcp45",
					"ssp126": "rcp26"}

		# Define the input file names based on scenario
		eais_filename = "dp21_eais_{}.nc".format(scen_dict[scenario])
		wais_filename = "dp21_wais_{}.nc".format(scen_dict[scenario])

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

		return years, eais_samps, wais_samps


def FindRefVals(timeseries, years, baseyear):

	# Append a zero to the beginning of the timeseries at year 2000
	# This was used for Bob's version of the DP20 data, not the current available data
	#timeseries = np.append(np.array([0.0]), timeseries)
	#years = np.append(np.array([2000]), years)

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
	parser = argparse.ArgumentParser(description="Run the DP21 icesheet pre-processing stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="Emissions scenario of interest [default = rcp85]", default="spp585")
	parser.add_argument('--baseyear', help="Year to which projections are referenced [default = 2000]", default=2000, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str, default="")

	# Parse the arguments
	args = parser.parse_args()

	# If baseyear is before year 2000, fail because we have no data prior to year 2000
	if args.baseyear < 2000:
		raise Exception("Base year must be greater than 2000: baseyear = {}".format(args.baseyear))

	# Run the preprocessing stage with the provided arguments
	dp21_preprocess_icesheet(args.scenario, args.baseyear, args.pipeline_id, args.climate_data_file)

	exit()
