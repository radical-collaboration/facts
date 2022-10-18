import os
import sys
import argparse
import pickle
import xarray as xr


def fair_fit_temperature(param_file, pipeline_id):

	# Load the AR6 calibrated parameters for the FAIR model
	pars = xr.load_dataset("./parameters/fair_ar6_climate_params_v4.0.nc")

	# Save the fit data to a pickle
	output = {"pars": pars, "param_file": param_file}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=-1)
	outfile.close()

	return(None)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fit stage for the FAIR AR6 temperature module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--param_file', help="Full path to FAIR parameter file", default="./parameters/fair_ar6_climate_params_v4.0.nc")

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	fair_fit_temperature(param_file=args.param_file, pipeline_id=args.pipeline_id)


	sys.exit()
