import os
import sys
import argparse
import pickle


def emulandice_fit_GrIS(pipeline_id):

	# Define the linear trend terms to be added to the samples
	trend_mean = {"GIS": 0.19}
	trend_sd = {"GIS": 0.1}

	# Populate the output dictionary
	outdata = {'trend_mean': trend_mean, 'trend_sd': trend_sd}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(outdata, outfile)
	outfile.close()

	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the fit stage for the emulandice GrIS module.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_fit_GrIS(args.pipeline_id)

	# Done
	sys.exit()
