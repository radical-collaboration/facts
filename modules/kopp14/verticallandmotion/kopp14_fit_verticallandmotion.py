import argparse

''' kopp14_fit_verticallandmotion.py

This runs the fitting stage for the vertical land motion component of the Kopp14
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Note: This is currently a NULL process. Results from the original Kopp14 workflow are
read in during the preprocessing. That data is then used directly to make projections in
the next stage. When time permits, we would like to have the full Kopp14 vertical land
motion component ported over.

'''

def kopp14_fit_verticallandmotion(pipeline_id):

	return(0)

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the Kopp14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_fit_verticallandmotion(args.pipeline_id)
	
	# Done
	exit()