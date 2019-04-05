import numpy as np
import pickle
import os
from IncludeModels import IncludeModels
from SmoothZOSTOGA import SmoothZOSTOGA
from DriftCorr import DriftCorr

''' kopp14_preprocess_thermalexp.py

This script runs the ocean dynamics pre-processing task for the Kopp 2014 workflow. 
This task generates the variables needed to configure the ocean dynamics submodel.

Parameters: 
None

Output: Pickle file containing the configuration for the ocean dynamics submodel

'''

def kopp14_preprocess_thermalexp(modeldir):
	
	# Define variables
	scenario = "rcp85"
	datayears = np.arange(1861,2300)
	targyears = np.arange(2010,2201,10)
	mergeZOSZOSTOGA = 1
	smoothwin = 19
	driftcorr = False
	baseyear = 2000
	GCMprobscale = 0.833
	
	# Read in the ZOSTOGA data
	modeldir = os.path.join(modeldir, scenario)
	(modellist, ZOSTOGA) = IncludeModels(modeldir, ("ZOSTOGA", "ZOSGA"), datayears)
	
	# Center, suture, and smooth ZOSTOGA
	sZOSTOGA = np.nan * ZOSTOGA
	for i in np.arange(0,ZOSTOGA.shape[1]):
		(ZOSTOGA[:,i], sZOSTOGA[:,i]) = SmoothZOSTOGA(ZOSTOGA[:,i], datayears, baseyear, smoothwin)
	
	# Apply the drift correction if needed
	if(driftcorr):
		(sZOSTOGA, CWdrift, histGICrate, selectyears) = DriftCorr(sZOSTOGA, datayears, baseyear, modellist)
	else:
		CWdrift = np.nan
		histGICrate = np.nan
		selectyears = np.nan
	
	# Store the configuration in a pickle
	output = {'scenario': scenario, 'datayears': datayears,\
		'targyears': targyears,	'mergeZOSZOSTOGA': mergeZOSZOSTOGA,\
		'smoothwin': smoothwin, 'driftcorr': driftcorr, 'baseyear': baseyear,\
		'GCMprobscale': GCMprobscale}
	
	# Write the configuration to a file
	outdir = os.path.join(os.path.dirname(__file__), "data")
	outfile = open(os.path.join(outdir, "kopp14_thermalexp_config.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Store the ZOSTOGA variables in a pickle
	output = {'ZOSTOGA': sZOSTOGA, 'modellist': modellist, 'CWdrift': CWdrift,\
		'histGICrate': histGICrate, 'selectyears': selectyears}
	
	# Write the ZOSTOGA variables to a file
	outfile = open(os.path.join(outdir, "kopp14_thermalexp_ZOSTOGA.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Pass the model directory in via command line		
	try:
		modeldir = sys.argv[1]
	except:
		modeldir = os.path.join(os.path.dirname(__file__), "data", "SLR_ALL")
	finally:
		kopp14_preprocess_thermalexp(modeldir)
	
	exit()