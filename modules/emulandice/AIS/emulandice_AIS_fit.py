import os
import sys
import pickle
import argparse


def emulandice_fit_AIS(pipeline_id):

	# Define the linear trend terms to be added to the samples
	trend_mean = {"EAIS": 0.09, "WAIS": 0.18, "PEN": 0.06}
	trend_sd = {"EAIS": 0.04, "WAIS": 0.09, "PEN": 0.03}

	# Populate the output dictionary
	outdata = {'trend_mean': trend_mean, 'trend_sd': trend_sd}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(outdata, outfile)
	outfile.close()

	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the fit stage for the emulandice AIS module.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_fit_AIS(args.pipeline_id)

	# Done
	sys.exit()
