import numpy as np
import pickle as p
import os
import sys
import argparse
import scipy.io

''' kopp14_preprocess_verticallandmotion.py

This runs the preprocessing stage for the vertical land motion component of the Kopp14
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Note: The script currently loads the results saved from the original ProjectSL code in
the Kopp14 workflow. These then get passed through a NULL process "fitting" stage and 
pick up at the global projection stage.  With time we would like to port the entirety of
the Kopp14 vertical land motion component.

'''

def kopp14_preprocess_verticallandmotion(pipeline_id):

	# Read in the matlab output file
	mat_file = os.path.join(os.path.dirname(__file__), "SLRProjections161027GRIDDEDcore.mat")
	mat = scipy.io.loadmat(mat_file)
	
	# Extract the relevant data
	targyears = mat['targyears'][0]
	targregions = mat['targregions']
	targregionnames = mat['targregionnames']
	rateprojs = mat['rateprojs']
	rateprojssd = mat['rateprojssd']
	
	# Define the base year
	baseyear = 2000
	
	# Populate the output dictionary
	outdata = {'targyears': targyears, 'targregions': targregions, \
		'targregionnames': targregionnames, 'rateprojs': rateprojs, \
		'rateprojssd': rateprojssd, 'baseyear': baseyear}
	
	
	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_preprocess_verticallandmotion(args.pipeline_id)
	
	# Done
	exit()