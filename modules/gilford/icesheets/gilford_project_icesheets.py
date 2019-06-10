import pickle
import os
import argparse
from sample_emulator import *

''' gilford_project_icesheets.py

This script runs the icesheet projection task for the Gilford icesheet emulator. 

Parameters: 
config_file = Pickled file with configuration information generated in pre-processing
n_ts = Number of timeseries samples
n_lhc = Number of latin hypercube samples
seed = Random number generator seed value

Output: Pickle file containing the projections for the ice sheet emulator

Note: The total number of samples will be (n_ts x n_lhc)

'''

def gilford_project_icesheets(config_file, n_ts, n_lhc, seed):
	
	# Load the configuration file
	try:
		f = open(config_file, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	train_data_path = my_config["train_data_path"]
	emulator_path = my_config["emulator_path"]
	
	# Generate the samples from the emulator
	(samps, years) = sample_emulator(train_data_path, emulator_path, n_lhc, n_ts, seed)
	
	# Store the configuration in a pickle
	output = {'samps': samps, 'years': years}
	
	# Write the configuration to a file
	outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
	outfile = open(os.path.join(outdir, "gilford_icesheets_projections.pkl"), 'wb')
	pickle.dump(output, outfile, protocol=2)
	#pickle.dump(output, outfile)
	outfile.close()


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the ice sheet projection stage for the Gilford Ice Sheet Emulator workflow",\
	epilog="Note: This is meant to be run as part of the Gilford Ice Sheet module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--n_ts', help="Number of samples across time [default=10]", default=10, type=int)
	parser.add_argument('--n_lhc', help="Number of latin hypercube samples [default=1000]", default=1000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	
	parser.add_argument('--config_file', help="Configuration file produced in the pre-processing stage",\
	default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "gilford_icesheets_config.pkl"))
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection stage
	gilford_project_icesheets(args.config_file, args.n_ts, args.n_lhc, args.seed)
	
	exit()