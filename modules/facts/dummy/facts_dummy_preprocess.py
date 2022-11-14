import argparse


''' facts_dummy_preprocess.py

This module does nothing, but provides a convenient module to use for
moving files around on the server as specified in a custom pipeline.yml
file. (For example, this could be used to place climate data files
from a previous run in the $SHARED/climate directory.)

'''

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the facts/dummy module stage.",\
	epilog="Note: This is meant to be run within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	exit()
