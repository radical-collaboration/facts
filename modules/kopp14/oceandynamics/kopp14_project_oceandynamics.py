import argparse

''' kopp14_project_oceandynamics.py

This runs the projection stage for the ocean dynamics component of the Kopp14
workflow. Projections for ocean dynamics are regional, so the global projection stage
is a NULL process.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code


'''

def kopp14_project_oceandynamics(pipeline_id):
	
	return(0)

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the Kopp14 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the project stage with the user defined RCP scenario
	kopp14_project_oceandynamics(args.pipeline_id)
	
	# Done
	exit()