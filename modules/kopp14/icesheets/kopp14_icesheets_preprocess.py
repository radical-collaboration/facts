import numpy as np
import pickle as p
import os
import sys
import argparse

''' kopp14_preprocess_icesheets.py

This script runs the ice sheet pre-processing task for the Kopp 2014 workflow. Currently,
the data are hard-coded in this script, from which a pickled file is generated. This
could easily be altered to take filenames as parameters from which these data are read.
This would be useful for testing multiple sets of these data (likely with different
assumptions) in batch.

Parameters:
rcp_scenario = RCP scenario (default = 'rcp85')
use_srocc = Use a calibration more closely in tune with SROCC
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickled file containing barates, lastdecadegt, and aris2090 variables within a
dictionary vairable.

'''

def kopp14_preprocess_icesheets(rcp_scenario, use_srocc, pipeline_id):

	# Fill the rates data matricies
	barates = np.array([[0.8, 1.0, 1.2, 2.4, 5.8],\
						[0.2, 3.0, 0.3, 1.5, 11.8],\
						[-1.9, 2.8, -1.5, 0.2, 10.2]])
	lastdecadegt = np.array([-211, -85-29, 26])

	# Standard K14 calibration
	aris2090 = np.array([[[0.07, 0.12, 0.21], [-0.06, 0.04, 0.12]],\
						[[0.04, 0.08, 0.13], [-0.04, 0.05, 0.13]],\
						[[0.04, 0.08, 0.13], [-0.04, 0.05, 0.13]],\
						[[0.04, 0.06, 0.10], [-0.03, 0.05, 0.14]]]) * 1000

	# Values consistent with SROCC 2019
	if use_srocc:
		aris2090 = np.array([[[0.07, 0.12, 0.21], [0.02, 0.10, 0.23]],\
						[[0.04, 0.08, 0.13], [0.01, 0.05, 0.13]],\
						[[0.04, 0.08, 0.13], [0.01, 0.05, 0.13]],\
						[[0.04, 0.06, 0.10], [0.01, 0.04, 0.10]]]) * 1000

	# Fill the correlation data matricies
	bacorris = np.array([[1.0, 0.7, -0.2], [0.7, 1.0, -0.2], [-0.2, -0.2, 1]])
	arcorris = np.array([[1.0, 0.5], [0.5, 1.0]])
	#bacorris = np.eye(3)
	#arcorris = np.eye(2)

	# Which RCP scenario does the user want?
	rcp_list = ["rcp85", "rcp60", "rcp45", "rcp26"]
	try:
		rcp_ind = rcp_list.index(rcp_scenario)
	except:
		print("kopp14 icesheets: Invalid RCP scenario \"{0}\"".format(rcp_scenario))
		sys.exit(1);

	# Collate rates data into a single dictionary
	data_is = {'barates': barates, 'lastdecadegt': lastdecadegt, 'aris2090': aris2090[rcp_ind,:,:]}

	# Collate the correlation data into a single dictionary
	corr_is = {'bacorris': bacorris, 'arcorris': arcorris, 'scenario': rcp_scenario}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	ratefile = open(os.path.join(outdir, "{}_rates.pkl".format(pipeline_id)), 'wb')
	p.dump(data_is, ratefile)
	ratefile.close()

	# Write the correlation data to a pickle file
	corrfile = open(os.path.join(outdir, "{}_corr.pkl".format(pipeline_id)), 'wb')
	p.dump(corr_is, corrfile)
	corrfile.close()


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', '-s', help="RCP Scenario", choices=['rcp85', 'rcp60', 'rcp45', 'rcp26'], default='rcp85')
	parser.add_argument('--srocc', help="Use the SROCC calibration (1=SROCC, 0=Standard; default=0)", choices=[0,1], type=int, default=0)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_preprocess_icesheets(args.scenario, args.srocc, args.pipeline_id)

	# Done
	exit()
