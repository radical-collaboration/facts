import numpy as np
import pickle
import os
import sys
import re
import argparse
from IncludeModels import IncludeModels
from IncludeDABZOSModels import *
from SmoothZOSTOGA import SmoothZOSTOGA
from DriftCorr import DriftCorr
from read_locationfile import ReadLocationFile
from Smooth import Smooth

''' kopp14_preprocess_oceandynamics.py

This runs the preprocessing stage for the ocean dynamics component of the Kopp14
workflow.

Parameters:
rcp_scenario = RCP scenario of interest (default = "rcp85")
zostoga_modeldir = Directory that contains the ZOSTOGA model data
zos_modeldir = Directory that contains the ZOS model data (in *.mat format)
driftcorr = Apply the drift correction?
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code


'''

def kopp14_preprocess_oceandynamics(rcp_scenario, zostoga_modeldir, zos_modeldir, driftcorr, baseyear, pyear_start, pyear_end, pyear_step, locationfilename, pipeline_id):

	# Define variables
	datayears = np.arange(1861,2300)
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	targyears = np.union1d(targyears, baseyear)
	mergeZOSZOSTOGA = True
	smoothwin = 19
	GCMprobscale = 0.833
	maxDOF = np.iinfo(np.int32).max

	#--------- Begin Thermal Expansion ---------------------------------------------------
	# Read in the ZOSTOGA data
	zostoga_modeldir = os.path.join(zostoga_modeldir, rcp_scenario)
	(zostoga_modellist, ZOSTOGA) = IncludeModels(zostoga_modeldir, ("ZOSTOGA", "ZOSGA"), datayears)

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
		'GCMprobscale': GCMprobscale, 'maxDOF': maxDOF}

	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_config.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Store the ZOSTOGA variables in a pickle
	output = {'sZOSTOGA': sZOSTOGA, 'zostoga_modellist': zostoga_modellist, 'CWdrift': CWdrift,\
		'histGICrate': histGICrate, 'selectyears': selectyears}

	# Write the ZOSTOGA variables to a file
	outfile = open(os.path.join(outdir, "{}_ZOSTOGA.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()



	#------------ Begin Ocean Dynamics ---------------------------------------------------

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, focus_site_ids, focus_site_lats, focus_site_lons) = ReadLocationFile(locationfile)

	# Load the ZOS data
	(zos_modellist, ZOS_raw) = IncludeDABZOSModels(zos_modeldir, rcp_scenario, focus_site_lats, focus_site_lons)

	# Find the overlap between ZOS and ZOSTOGA
	comb_modellist, zostoga_model_idx, zos_model_idx = np.intersect1d(zostoga_modellist, zos_modellist, return_indices=True)

	'''
	NOTE: POTENTIAL BUG IN ORIGINAL CODE
	The original code uses the smoothed ZOSTOGA data as the raw ZOSTOGA values. This in
	turn smooths the ZOSTOGA data over the 'smoothwin' period twice. To replicate the
	bug, set 'ZOSTOGAadj' to a subset of 'sZOSTOGA' instead of 'ZOSTOGA'.
	'''
	ZOSTOGAadj = sZOSTOGA[:,zostoga_model_idx]  # Replicate potential bug
	#ZOSTOGAadj = ZOSTOGA[:,zostoga_model_idx]  # Fix for potential bug
	ZOS_raw = ZOS_raw[:, zos_model_idx, :]

	# Should we merge ZOSTOGA and ZOS?
	# NOTE: ZOS starts at 1860 while ZOSTOGA starts at 1861. Pop off the 1860 value for
	# ZOS and merge if necessary
	if(mergeZOSZOSTOGA):
		ZOS = ZOS_raw + ZOSTOGAadj[:,:,np.newaxis]
	else:
		ZOS = ZOS_raw

	# Smooth ZOS and ZOSTOGA over 19 year smoothing window
	def nanSmooth(x, w):
		idx = np.flatnonzero(~np.isnan(x))
		temp = x
		if len(idx) > 1:
			temp[idx] = Smooth(x[idx], w)
		return(temp)

	sZOS = np.apply_along_axis(nanSmooth, axis=0, arr=ZOS, w=smoothwin)
	sZOSTOGAadj = np.apply_along_axis(nanSmooth, axis=0, arr=ZOSTOGAadj, w=smoothwin)

	# Center the smoothed ZOS/ZOSTOGAadj to the baseyear
	baseyear_idx = np.flatnonzero(datayears == baseyear)
	sZOS = np.apply_along_axis(lambda z, idx: z - z[idx], axis=0, arr=sZOS, idx=baseyear_idx)
	sZOSTOGAadj = np.apply_along_axis(lambda z, idx: z - z[idx], axis=0, arr=sZOSTOGAadj, idx=baseyear_idx)

	# Store the ZOS variable in a pickle
	output = {'sZOS': sZOS, 'zos_modellist': zos_modellist, 'datayears': datayears, \
				'focus_site_ids': focus_site_ids, 'focus_site_lats': focus_site_lats, \
				'focus_site_lons': focus_site_lons, 'sZOSTOGAadj': sZOSTOGAadj, 'comb_modellist': comb_modellist}

	# Write the ZOS variables to a file
	outfile = open(os.path.join(outdir, "{}_ZOS.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp14 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="RCP Scenario [default=\'rcp85\']", choices=['rcp85', 'rcp60', 'rcp45', 'rcp26'], default='rcp85')

	parser.add_argument('--zostoga_model_dir', help="Directory containing ZOSTOGA GCM output",\
	default=os.path.join(os.path.dirname(__file__), "SLR_ALL"))

	parser.add_argument('--zos_model_dir', help="Directory containing ZOS GCM output",\
	default=os.path.join(os.path.dirname(__file__), "zosfinal"))

	parser.add_argument('--no_drift_corr', help="Do not apply the drift correction", action='store_true')

	parser.add_argument('--baseyear', help="Base or reference year for projetions [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)

	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")

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


	# Pass the model directory in via command line
	kopp14_preprocess_oceandynamics(args.scenario, args.zostoga_model_dir, args.zos_model_dir, not args.no_drift_corr, args.baseyear, args.pyear_start, args.pyear_end, args.pyear_step, args.locationfile, args.pipeline_id)

	# Done
	exit()
