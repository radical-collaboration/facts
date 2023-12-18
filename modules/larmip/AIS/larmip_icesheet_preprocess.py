import numpy as np
import argparse
import h5py
import os
import sys
import pickle


def GetSATData(fname, scenario, refyear_start=1850, refyear_end=1900, year_start=1900, year_end=2300):

	# Open the file
	df_ssp = h5py.File(fname,'r')

	# Extract the surface temperature for this scenario
	if scenario not in df_ssp.keys():
		raise Except("Scenario {} not found in {}".format(scenario, fname))
	ssp_folder = df_ssp.get(scenario)
	sat_ssp = ssp_folder.get('surface_temperature')

	# Get the number of ensemble members from the data
	_,nens = sat_ssp.shape

	# Extract the years available
	sat_years = df_ssp['year'][()]

	# Which indices align with the reference and trim years
	refyear_start_idx = np.flatnonzero(sat_years == refyear_start)[0]
	refyear_end_idx = np.flatnonzero(sat_years == refyear_end)[0] #+ 1
	year_start_idx = np.flatnonzero(sat_years == year_start)[0]
	year_end_idx = np.flatnonzero(sat_years == year_end)[0] + 1

	# Normalize and crop temperature series
	Time = np.arange(year_start, year_end+1)
	SATave = np.mean(sat_ssp[refyear_start_idx:refyear_end_idx,:], axis=0)
	SAT = sat_ssp[year_start_idx:year_end_idx,:] - SATave

	# Close the h5 file
	df_ssp.close()

	# Done
	return(SAT, Time, nens)


def larmip_preprocess_icesheet(scenario, pipeline_id, fname):

	# Load the temperature data
	SAT,Time,NumTensemble = GetSATData(fname, scenario)
	#SAT,Time,NumTensemble = GetSATData("./data/IPCC_SSP_Projections/twolayer_SSPs.h5", scenario)
	Tlen = len(Time)

	# Save the preprocessed data to a pickle
	output = {"SAT": SAT, "Time": Time, "NumTensemble": NumTensemble, "Tlen": Tlen, "scenario": scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_preprocess.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(None)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the preprocessing stage for the LARMIP icesheet module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--scenario', help="SSP Emissions scenario", required=True)
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)

	# Parse the arguments
	args = parser.parse_args()

	# Run the larmip code
	larmip_preprocess_icesheet(pipeline_id=args.pipeline_id, scenario=args.scenario, fname=args.climate_data_file)


	sys.exit()
