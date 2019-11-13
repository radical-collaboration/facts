import numpy as np
import pickle
import os
import argparse
from IncludeModels import IncludeModels
from SmoothZOSTOGA import SmoothZOSTOGA
from DriftCorr import DriftCorr

''' kopp14_preprocess_thermalexp.py

This script runs the ocean dynamics pre-processing task for the Kopp 2014 workflow. 
This task generates the variables needed to configure the ocean dynamics submodel.

Parameters: 
rcp_scenario = RCP scenario of interest (default = "rcp85")
modeldir = Directory that contains the GCM model data
driftcorr = Apply the drift correction?
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickle file containing the configuration for the ocean dynamics submodel

'''

def kopp14_preprocess_thermalexp(rcp_scenario, modeldir, driftcorr, pipeline_id):
	
	# Define variables
	datayears = np.arange(1861,2300)
	targyears = np.arange(2010,2201,10)
	mergeZOSZOSTOGA = 1
	smoothwin = 19
	baseyear = 2005
	GCMprobscale = 0.833
	
	# Read in the ZOSTOGA data
	modeldir = os.path.join(modeldir, rcp_scenario)
	(modellist, ZOSTOGA) = IncludeModels(modeldir, ("ZOSTOGA", "ZOSGA"), datayears)
	
	# Center, suture, and smooth ZOSTOGA
	sZOSTOGA = np.nan * ZOSTOGA
	for i in np.arange(0,ZOSTOGA.shape[1]):
		(ZOSTOGA[:,i], sZOSTOGA[:,i]) = SmoothZOSTOGA(ZOSTOGA[:,i], datayears, baseyear, smoothwin)
	
	# Apply the drift correction if needed
	if(driftcorr):
		gslfile = os.path.join(os.path.dirname(__file__), "CSIRO_Recons_gmsl_yr_2011.csv")
		(sZOSTOGA, CWdrift, histGICrate, selectyears) = DriftCorr(sZOSTOGA, datayears, baseyear, rcp_scenario, gslfile)
	else:
		CWdrift = np.nan
		histGICrate = np.nan
		selectyears = np.nan
	
	# Store the configuration in a pickle
	output = {'rcp_scenario': rcp_scenario, 'datayears': datayears,\
		'targyears': targyears,	'mergeZOSZOSTOGA': mergeZOSZOSTOGA,\
		'smoothwin': smoothwin, 'driftcorr': driftcorr, 'baseyear': baseyear,\
		'GCMprobscale': GCMprobscale}
	
	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_config.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Store the ZOSTOGA variables in a pickle
	output = {'ZOSTOGA': sZOSTOGA, 'modellist': modellist, 'CWdrift': CWdrift,\
		'histGICrate': histGICrate, 'selectyears': selectyears}
	
	# Write the ZOSTOGA variables to a file
	outfile = open(os.path.join(outdir, "{}_ZOSTOGA.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion pre-processing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="RCP Scenario [default=\'rcp85\']", choices=['rcp85', 'rcp60', 'rcp45', 'rcp26'], default='rcp85')
	
	parser.add_argument('--model_dir', help="Directory containing ZOS/ZOSTOGA GCM output",\
	default=os.path.join(os.path.dirname(__file__), "SLR_ALL"))
	
	parser.add_argument('--no_drift_corr', help="Do not apply the drift correction", action='store_true')
	
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Pass the model directory in via command line
	kopp14_preprocess_thermalexp(args.scenario, args.model_dir, not args.no_drift_corr, args.pipeline_id)
	
	exit()