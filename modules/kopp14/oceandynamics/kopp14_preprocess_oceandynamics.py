import numpy as np
import pickle as p
import os
import sys
import argparse
import scipy.io

''' kopp14_preprocess_oceandynamics.py

This runs the preprocessing stage for the ocean dynamics component of the Kopp14
workflow.

Parameters:
rcp_scenario = RCP scenario (default = 'rcp85')
pipeline_id = Unique identifier for the pipeline running this code

Note: The script currently loads the results saved from the original ProjectSL code in
the Kopp14 workflow. These then get passed through a NULL process "fitting" stage and 
pick up at the global projection stage.  With time we would like to port the entirety of
the Kopp14 ocean dynamics component.

'''

def kopp14_preprocess_oceandynamics(rcp_scenario, pipeline_id):

	# Which RCP scenario does the user want?
	rcp_list = ["rcp85", "rcp60", "rcp45", "rcp26"]
	try:
		rcp_ind = rcp_list.index(rcp_scenario)
	except:
		print("kopp14 oceandynamics: Invalid RCP scenario \"{0}\"".format(rcp_scenario))
		sys.exit(1);

	# Read in the matlab output file
	mat_file = os.path.join(os.path.dirname(__file__), "SLRProjections161027GRIDDEDcore.mat")
	mat = scipy.io.loadmat(mat_file)
	
	# Extract the relevant data
	targyears = np.squeeze(mat['targyears'])
	targregions = mat['targregions']
	targregionnames = mat['targregionnames']
	OceanDynYears = mat['OceanDynYears']
	ThermExpYears = mat['ThermExpYears']
	OceanDynMean = np.squeeze(mat['OceanDynMean'][:,:,rcp_ind])
	OceanDynStd = np.squeeze(mat['OceanDynStd'][:,:,rcp_ind])
	OceanDynN = np.squeeze(mat['OceanDynN'][:,:,rcp_ind])
	OceanDynTECorr = np.squeeze(mat['OceanDynTECorr'][:,:,rcp_ind])
	tesamps = np.squeeze(mat['samps'][:,mat['colTE']-1,:,rcp_ind])
	ThermExpMean = np.squeeze(mat['ThermExpMean'][:,rcp_ind])
	ThermExpStd = np.squeeze(mat['ThermExpStd'][:,rcp_ind])
	
	# Define the GCM probability scale
	GCMprobscale = 0.833
	
	# Define the maximum degrees of freedom
	maxDOF = np.inf
	
	# Populate the output dictionary
	outdata = {'targyears': targyears, 'OceanDynYears': OceanDynYears, 'OceanDynN': OceanDynN, \
		'ThermExpYears': ThermExpYears, 'OceanDynMean': OceanDynMean, 'OceanDynStd': OceanDynStd, \
		'OceanDynTECorr': OceanDynTECorr, 'tesamps': tesamps, 'ThermExpMean': ThermExpMean, \
		'ThermExpStd': ThermExpStd, 'targregions': targregions, 'targregionnames': targregionnames, \
		'GCMprobscale': GCMprobscale, 'maxDOF': maxDOF}
	
	
	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp14 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', '-s', help="RCP Scenario", choices=['rcp85', 'rcp60', 'rcp45', 'rcp26'], default='rcp85')
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_preprocess_oceandynamics(args.scenario, args.pipeline_id)
	
	# Done
	exit()