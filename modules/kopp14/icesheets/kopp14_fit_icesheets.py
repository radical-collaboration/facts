import pickle
import sys
import os
import argparse
import numpy as np
from CalcISDists import CalcISDists as calcDists

''' kopp14_fit_icesheets.py

This script runs the ice sheet fitting task for the Kopp 2014 workflow. Data produced
from the pre-processing task 'kopp14_preprocess_icesheets.py' are used here. This task 
requires a single paramter, which is the name of the pickled file produced by the 
pre-processing task. Output is another pickled file that contains the fitted parameters
for the log-normal distributions representing ice sheet contributions to future sea-level 
rise

Parameters: 
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickled file containing the fitted parameter values


'''

def kopp14_fit_icesheets(pipeline_id):
	
	# Open the file
	infile = "{}_rates.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
	
	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()
	
	# Fit the distributions
	(batheteais, bathetwais, bathetgis, arthetais, arthetgis, islastdecade) = calcDists(my_data['barates'], my_data['lastdecadegt'], my_data['aris2090'])
	
	# Put the results into a dictionary
	output = {'batheteais': batheteais, 'bathetwais': bathetwais, 'bathetgis': bathetgis,\
				'arthetais': arthetais, 'arthetgis': arthetgis,\
				'islastdecade': islastdecade}
	
	# Write the results to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the fitting process using the provided rate file from preprocessing stage
	kopp14_fit_icesheets(args.pipeline_id)
	
	# Done
	exit()