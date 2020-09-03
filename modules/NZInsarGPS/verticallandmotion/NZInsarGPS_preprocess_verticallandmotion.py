import numpy as np
import pickle as p
import os
import sys
import argparse

''' NZInsarGPS_preprocess_verticallandmotion.py

This runs the preprocessing stage for the vertical land motion component of the NZInsarGPS
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code
baseyear = The year from which projections should be zeroed

'''

def DataCheck(x):
	if("NaN" in x.strip()):
		return(np.nan)
	else:
		return(float(x))


def ReadData(filename):
	
	# Initialize variables to hold data and site information
	lats = []
	lons = []
	rates = []
	sds = []
	ids = []
	
	# Initialize IDs
	this_id = 0
	
	# Define the rate file name
	rate_file = os.path.join(os.path.dirname(__file__), filename)
	
	# Open the rate file
	with open(rate_file, 'r') as f:
		
		# Skip the first line (header info)
		line = f.readline()
		
		# Loop over the lines of the file
		for line in f:
			
			# Increment the ID
			this_id = this_id + 1
			
			# Split the line into components
			line_parts = [str(x) for x in line.split("\t")]
			(this_lon, this_lat, this_rate, this_sd) = line_parts
			
			# Test if any data is available for this entry
			if np.isnan(DataCheck(this_rate)):
				continue
			
			# Store the information
			ids.append(this_id)
			lats.append(DataCheck(this_lat))
			lons.append(DataCheck(this_lon))
			rates.append(DataCheck(this_rate))
			sds.append(DataCheck(this_sd))
			
	# Convert lists into numpy arrays and remove entries with no data
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	rates = np.array(rates)
	sds = np.array(sds)
	
	return(ids, lats, lons, rates, sds)


def NZInsarGPS_preprocess_verticallandmotion(pipeline_id, baseyear, inputtype):
	
	# Define the input file based on input type
	type2file = {"grid_insar": "VLM-grid_insar.dat", "grid_gps": "VLM-grid_gps.dat", "hires": "coast_ins.dat"}
	if inputtype not in type2file.keys():
		raise Exception("Input type not recognized: {}".format(inputtype))
	inputfile = type2file[inputtype]

	# Read input file
	(ids, lats, lons, rates, sds) = ReadData(inputfile)
	
	# Populate the output dictionary
	outdata = {'ids': ids, 'lats': lats, 'lons': lons, 'rates': rates, 'sds': sds,\
				'baseyear': baseyear, 'nsites': len(ids), 'inputfile': inputfile}
	
	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the NZInsarGPS vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the NZInsarGPS module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--baseyear', help="Baseline year from which to calculate changes in sea level", default=2005, type=int)
	parser.add_argument('--inputtype', help="Type of data to use (default=\"grid_insar\")", choices=["grid_insar", "grid_gps", "hires"], default="grid_insar")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	NZInsarGPS_preprocess_verticallandmotion(args.pipeline_id, args.baseyear, args.inputtype)
	
	# Done
	exit()