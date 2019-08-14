import numpy as np
import pickle
import os
import argparse

''' kopp14_fit_glaciers.py

This script runs the glaciers fitting task for the Kopp 2014 workflow. 
This task generates the t-distributions for projecting glacier contributions.

Parameters: 
data_file = The glacier data/configuration file

Output: 

'''

def kopp14_fit_glaciers(data_file):
	
	# Load the configuration file
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	rcp_scenario = my_config["rcp_scenario"]
	projGIC = my_config["projGIC"]
	projGICse = my_config["projGICse"]
	projGICyrs = my_config["projGICyrs"]
	projGICmodel = my_config["projGICmodel"]
	
	
	## Store the thermal expansion variables in a pickle
	#output = {'ThermExpMean': ThermExpMean, 'ThermExpStd': ThermExpStd,\
	#	'ThermExpYears': ThermExpYears,	'ThermExpN': ThermExpN}
	#outfile = open(os.path.join(os.path.dirname(__file__), "kopp14_thermalexp_fit.pkl"), 'wb')
	#pickle.dump(output, outfile)
	#outfile.close()


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glaciers fitting stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--data_file', help="Pickle file containing Glacier data produced from preprocessing stage [default=\'kopp14_glaciers_data.pkl\']", default='kopp14_glaciers_data.pkl')
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the fitting process with the provided command line arguments
	kopp14_fit_glaciers(args.data_file)
	
	exit()