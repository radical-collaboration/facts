import numpy as np
import pickle as p
import os
import sys
import argparse
from netCDF4 import Dataset

''' caron18_GIA_preprocess.py

This runs the preprocessing stage for the GIA contribution to relative sea-level change from Caron et al. (2018).


Parameters:
pipeline_id = Unique identifier for the pipeline running this code
baseyear = The year from which projections should be zeroed

'''

def caron18_preprocess_GIA(pipeline_id):

	# Initialize variables to hold data and site information
	lats = []
	lons = []
	rates = []
	sds = []

	# Define the rate file name
	rate_file = os.path.join(os.path.dirname(__file__), "GIA_stats.nc")
	ds = Dataset(rate_file)

	# Extract the relevant data
	lat = np.array(ds['lat'][:])
	lon = np.array(ds['lon'][:])

	lons, lats = np.meshgrid(lon, lat)
	lats = lats.flatten()
	lons = lons.flatten()

	rates = np.array(ds['rsl_mean'][:]).flatten()
	sds = np.array(ds['rsl_sterr'][:]).flatten()

	# Populate the output dictionary
	outdata = { 'lats': lats, 'lons': lons, 'rates': rates,\
				'sds': sds}


	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Caron 18 GIA workflow",\
	epilog="Note: This is meant to be run as part of the GIA module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	caron18_preprocess_GIA(args.pipeline_id)

	# Done
	exit()
