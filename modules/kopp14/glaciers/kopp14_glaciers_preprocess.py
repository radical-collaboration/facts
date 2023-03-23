import numpy as np
import pickle
import os
import argparse
from readMarzeion import readMarzeion

''' kopp14_preprocess_glaciers.py

This script runs the glaciers pre-processing task for the Kopp 2014 workflow. 
This task generates the variables needed to configure the glaciers submodel.

Parameters: 
rcp_scenario = The RCP scenario
baseyear = The reference year to which projections are zeroed
pyear_start = Year at which projections should begin
pyear_end = Year at which projections should end
pyear_start = Step size in years between pyear_start and pyear_end
pipeline_id = Unique identifier for the pipeline running this code

Output:
"%PIPELINE_ID%_data.pkl" = Contains the GIC data
"%PIPELINE_ID%_config.pkl" = Contains the configuration parameters
"%PIPELINE_ID%_fp.pkl" = Contains the fingerprint information

'''

def kopp14_preprocess_glaciers(rcp_scenario, baseyear, pyear_start, pyear_end, pyear_step, pipeline_id):
	
	# Use readMarzeion to read in the glacier data
	glacdir = "."
	fpmap = os.path.join(os.path.dirname(__file__), "fingerprint_region_map.csv")
	(projGIC, projGICse, projGICyrs, projGICmodel, fpmapperids, fpmaps, _) = readMarzeion(rcp_scenario, glacdir, fpmap, baseyear, discardAntarctica=True)
	
	# Define the target years for projections
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	
	# Insert the base year into the target years
	targyears = np.union1d(targyears, baseyear)
	
	# Store the data in a pickle
	output = {'projGIC': projGIC, 'projGICse': projGICse, 'projGICyrs': projGICyrs, 'projGICmodel': projGICmodel}
	
	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Store the configuration in a pickle
	output = {'rcp_scenario': rcp_scenario, 'baseyear': baseyear, 'targyears': targyears}
	
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
	parser.add_argument('--baseyear', help="Base or reference year for projetions [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Make sure the base year and target years are within data limits for this module
	if(args.baseyear < 2000):
		raise Exception("Base year cannot be less than year 2000: baseyear = {}".format(args.baseyear))
	if(args.baseyear > 2300):
		raise Exception("Base year cannot be greater than year 2300: baseyear = {}".format(args.baseyear))
	if(args.pyear_start < 2000):
		raise Exception("Projection year cannot be less than year 2000: pyear_start = {}".format(args.pyear_start))
	if(args.pyear_end > 2300):
		raise Exception("Projection year cannot be greater than year 2300: pyear_end = {}".format(args.pyear_end))
	
	# Make sure the target year stepping is positive
	if(args.pyear_step < 1):
		raise Exception("Projection year step must be greater than 0: pyear_step = {}".format(args.pyear_step))
	
	# Run the preprocessing stage with the provided arguments
	kopp14_preprocess_glaciers(args.scenario, args.baseyear, args.pyear_start, args.pyear_end, args.pyear_step, args.pipeline_id)
	
	exit()