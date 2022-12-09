import argparse
import sys

''' dp_fit_icesheet.py

Fitting process for the dp16 icesheet module.

There's no fitting process for the DP16 module since it's sampling from discrete
samples.  This is essentially a NULL process.

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def dp16_fit_icesheet(pipeline_id):
	
	return(None)


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the DP16 icesheet fitting stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	dp16_fit_icesheet(args.pipeline_id)
	
	sys.exit()