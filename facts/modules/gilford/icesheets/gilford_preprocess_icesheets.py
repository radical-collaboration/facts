import numpy as np
import pickle
import os
import argparse

''' gilford_preprocess_icesheets.py

This script runs the icesheet pre-processing task for the Gilford icesheet emulator. 

Parameters: 
train_data_path = Path to the training data
emulator_path = Path to the pre-generated emulator object

Output: Pickle file containing the configuration for the ice sheet emulator

'''

def gilford_preprocess_icesheets(train_data_path, emulator_path):
	
	# Define variables
	
	
	# Store the configuration in a pickle
	output = {'train_data_path': train_data_path, 'emulator_path': emulator_path}
	
	# Write the configuration to a file
	outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
	outfile = open(os.path.join(outdir, "gilford_icesheets_config.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the ice sheet pre-processing stage for the Gilford Ice Sheet Emulator workflow",\
	epilog="Note: This is meant to be run as part of the Gilford Ice Sheet module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--train_data_path', help="Path to the training data",\
	default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "d19_rcp_model_2100_traindata.pk1"))
	
	parser.add_argument('--emulator_path', help="Path to the pre-generated emulator objects",\
	default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "d19_rcp_model_2100.ckpt"))
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the pre-processing stage for the Gilford emulator workflow
	gilford_preprocess_icesheets(args.train_data_path, args.emulator_path)
	
	exit()