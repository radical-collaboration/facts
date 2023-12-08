import numpy as np
import os
import sys
import pickle
import argparse
import re

def tlm_fit_thermalexpansion(pipeline_id):

	# Load the preprocessed data
	data_file = "{}_tlmdata.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")

	# Extract the data variables
	my_data = pickle.load(f)
	f.close()

	expcoefs = my_data['expcoefs']
	rmses = my_data['rmses']
	expcoefs_models = my_data['expcoefs_models']
	rmses_models = my_data['rmses_models']

	# Cumulative probability at which to trim the distribution of models
	pb_clip = 0.85

	# Clip gsat rmse cdf
	rmses_sorted = np.sort(rmses)
	sort_idx = np.argsort(rmses)
	pbs = np.arange(0,rmses_sorted.size)/rmses_sorted.size

	# Find which models to cut out
	include_models = rmses_models[sort_idx][pbs <= pb_clip]

	# Extract the coefficients for the models that made the cut
	model_idx = np.isin(expcoefs_models, include_models)
	expcoefs_clipped = expcoefs[model_idx]

	# Extract normal distribution modes
	mean_expcoefs = np.mean(expcoefs_clipped)
	std_expcoefs = np.std(expcoefs_clipped)

	# Store preprocessed data in pickles
	output = {'mean_expcoefs': mean_expcoefs, 'std_expcoefs': std_expcoefs, \
				'include_models': include_models, 'pb_clip': pb_clip}

	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_tlmfit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile,protocol=4)
	outfile.close()

	return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the TLM ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Define the names of the intermediate files this run will generate
	outdir = os.path.dirname(__file__)
	tlmfitfile = os.path.join(outdir, "{}_tlmfit.pkl".format(args.pipeline_id))


	# Runs the TLM fitting stage if intermediate files are not present
	if os.path.isfile(tlmfitfile):
		print("{} found, skipping TE fitting".format(tlmfitfile))
	else:
		tlm_fit_thermalexpansion(args.pipeline_id)






	# Done
	sys.exit()
