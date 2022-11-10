import argparse

''' bamber19_fit_icesheets.py

This runs the fitting stage for the Bamber et al 2019 ice sheet component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Note: This is currently a NULL process.

'''

def bamber19_fit_icesheets(pipeline_id):

	return(0)

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the IPCC AR6 Bamber et al. 2019 ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the ipccar6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	bamber19_fit_icesheets(args.pipeline_id)
	
	# Done
	exit()