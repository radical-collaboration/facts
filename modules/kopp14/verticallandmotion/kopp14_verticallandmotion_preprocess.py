import numpy as np
import pickle as p
import os
import sys
import argparse

''' kopp14_preprocess_verticallandmotion.py

This runs the preprocessing stage for the vertical land motion component of the Kopp 14
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code
baseyear = The year from which projections should be zeroed

'''

def kopp14_preprocess_verticallandmotion(pipeline_id):

	# Initialize variables to hold data and site information
	names = []
	ids = []
	lats = []
	lons = []
	rates = []
	sds = []

	# Define the rate file name
	rate_file = os.path.join(os.path.dirname(__file__), "bkgdrate-210306.tsv")

	# Open the rate file
	with open(rate_file, 'r') as f:

		# Skip the first line (header info)
		line = f.readline()

		# Loop over the lines of the file
		for line in f:

			# Split the line into components
			(this_name, this_id, this_lat, this_lon, this_rate, this_sd) = line.split("\t")

			# Store the information
			names.append(this_name)
			ids.append(int(this_id))
			lats.append(float(this_lat))
			lons.append(float(this_lon))
			rates.append(float(this_rate))
			sds.append(float(this_sd))

	# Populate the output dictionary
	outdata = {'names': names, 'ids': ids, 'lats': lats, 'lons': lons, 'rates': rates,\
				'sds': sds}


	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp 14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp 14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_preprocess_verticallandmotion(args.pipeline_id)

	# Done
	exit()
