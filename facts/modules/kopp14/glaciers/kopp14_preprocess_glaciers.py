import numpy as np
import pickle
import os
import argparse
from readMarzeion import readMarzeion

''' kopp14_preprocess_glaciers.py

This script runs the glaciers pre-processing task for the Kopp 2014 workflow. 
This task generates the variables needed to configure the glaciers submodel.

Parameters: 
rcp_scenario = The RCP scenario (default: rcp85)
pipeline_id = Unique identifier for the pipeline running this code

Output:
"%PIPELINE_ID%_data.pkl" = Contains the GIC data
"%PIPELINE_ID%_config.pkl" = Contains the configuration parameters
"%PIPELINE_ID%_fp.pkl" = Contains the fingerprint information

'''

def kopp14_preprocess_glaciers(rcp_scenario, pipeline_id):
	
	# Use readMarzeion to read in the glacier data
	glacdir = "."
	fpmap = os.path.join(os.path.dirname(__file__), "fingerprint_region_map.csv")
	(projGIC, projGICse, projGICyrs, projGICmodel, fpmapperids, fpmaps, _) = readMarzeion(rcp_scenario, glacdir, fpmap, discardAntarctica=True)
	
	# Define the target years for projections
	targyears = np.arange(2010, 2101, 10)
	
	# Store the data in a pickle
	output = {'projGIC': projGIC, 'projGICse': projGICse, 'projGICyrs': projGICyrs, 'projGICmodel': projGICmodel}
	
	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Store the configuration in a pickle
	output = {'rcp_scenario': rcp_scenario, 'targyears': targyears}
	
	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_config.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Store the fingerprint variables in a pickle
	output = {'fpmapperids': fpmapperids, 'fpmaps': fpmaps}
	
	# Write the fingerprint variables to a file
	outfile = open(os.path.join(outdir, "{}_fp.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glaciers pre-processing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="RCP Scenario [default=\'rcp85\']", choices=['rcp85', 'rcp60', 'rcp45', 'rcp26'], default='rcp85')
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	kopp14_preprocess_glaciers(args.scenario, args.pipeline_id)
	
	exit()