import argparse

''' caron18_GIA_project.py

This runs the projection stage for the GIA contribution to relative sea-level change from Caron et al. (2018).


Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Note: This is currently a NULL process. Rates are already read in during the preprocess
stage.  Projections for individual locations and grid points are handled in
post-processing.

'''

def caron18_project_GIA(pipeline_id):

	return(0)

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the Caron 18 GIA workflow",\
	epilog="Note: This is meant to be run as part of the Caron 18 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	caron18_project_GIA(args.pipeline_id)
	
	# Done
	exit()