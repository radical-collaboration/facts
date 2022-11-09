import numpy as np
import os
import sys
import fnmatch
import argparse
import re
import pickle
import shutil
from netCDF4 import Dataset



def GetSamples(ncfile, years, baseyear):

	# Load the nc file
	with Dataset(ncfile, "r") as nc:

		# Extract the variables
		ncyears = nc.variables["years"][...]
		samples = np.squeeze(nc.variables["surface_temperature"][...])
		scenario = nc.getncattr("Scenario")

	# Extract the reference value
	baseyear_idx = np.flatnonzero(ncyears == baseyear)
	ref_val = samples[:,baseyear_idx]

	# Indices of years to extract and return
	year_idx = np.flatnonzero(np.isin(ncyears, years))

	# Squeeze out the location dimension (should be global temperature trajectories)
	samples = samples[:,year_idx] - ref_val

	# Return the samples
	return(samples, scenario)


def WriteToCSV(outfile, samples, mode="w"):

	# Open the csv file
	with open(outfile, mode) as f:

		# Loop through the samples
		for i in np.arange(samples.shape[0]):

			# Put together the output string
			out_string = ",".join([str(x) for x in ["FAIR", "FAIR_{}".format(i+1), "FACTS", *samples[i,:]]])

			# Write this string to the output file
			f.write(out_string)
			f.write("\n")

	# Done
	return(None)


def emulandice_preprocess(infile, baseyear, pipeline_id):

	# If no input file was passed, look for one produced by a pre-projection workflow
	if infile is None:
		indir = os.path.dirname(__file__)
		if indir == "":
			indir = "."
		files = fnmatch.filter(os.listdir(indir), "*_gsat.nc")
		if files is None:
			raise Exception("Unable to find a usable input file")
		infile = files[0]

	# Years
	years = np.arange(2015,2101)

	# Get the samples
	samps, scenario = GetSamples(infile, years, baseyear)

	# How many samples are we running?
	nsamps = samps.shape[0]

	# Append these samples to the output file
	headfile =  os.path.join(os.path.dirname(__file__), "FACTS_CLIMATE_FORCING.csv.head")
	outfile = os.path.join(os.path.dirname(__file__), "FACTS_CLIMATE_FORCING.csv")
	shutil.copyfile(headfile,outfile)
	WriteToCSV(outfile, samps, mode="a")

	# Save the preprocessed data to a pickle
	output = {"scenario": scenario, "baseyear": baseyear, "infile": infile, \
			"facts_data_file": outfile, "nsamps": nsamps}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_preprocess.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=-1)
	outfile.close()


	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the preprocess stage for the emulandice module.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--input_data_file', help="Full path for temperature trajectory input file", default=None)
	parser.add_argument('--baseyear', help="Year to which projections should be referenced", type=int, default=2005)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_preprocess(args.input_data_file, args.baseyear, args.pipeline_id)

	# Done
	sys.exit()
