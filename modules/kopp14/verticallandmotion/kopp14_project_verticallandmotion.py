import argparse

''' ipccar6_project_verticallandmotion.py

This runs the projection stage for the vertical land motion component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Note: This is currently a NULL process. Rates are already read in during the preprocess
stage.  Projections for individual locations and grid points are handled in
post-processing.

'''

def ipccar6_project_verticallandmotion(pipeline_id):

	return(0)

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the IPCC AR6 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the IPCC AR6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	ipccar6_project_verticallandmotion(args.pipeline_id)
	
	# Done
	exit()