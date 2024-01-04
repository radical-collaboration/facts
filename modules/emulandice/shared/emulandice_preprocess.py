import numpy as np
import argparse
import h5py
import os
import sys
import pickle
import shutil


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


def WriteToCSV(outfile, samples, mode="w"):

	# Open the csv file
	with open(outfile, mode) as f:

		# Loop through the samples
		for i in np.arange(samples.shape[0]):

			# Put together the output string
			out_string = ",".join([str(x) for x in ["FAIR", "FAIR_{}".format(i+1), "FACTS", *samples[i,:]]])

			# Write this string to the output file
			f.write(out_string)
			f.write("\n")

	# Done
	return(None)


def emulandice_preprocess(pipeline_id, scenario, fname, baseyear, pipeline_id):

	# Years
	years = np.arange(2015,2101)

	# Load the temperature data
	samps,Time,nsamps = GetSATData(fname, scenario)
	Tlen = len(Time)

	# emulandice-specific adjustment

	# Extract the reference value
	baseyear_idx = np.flatnonzero(Time == baseyear)
	ref_val = samps[baseyear_idx,:]

	# Indices of years to extract and return
	year_idx = np.flatnonzero(np.isin(Time, years))

	# Squeeze out the location dimension (should be global temperature trajectories)
	samps = np.transpose(samps[year_idx,:] - ref_val)

	# Append these samples to the output file
	headfile =  os.path.join(os.path.dirname(__file__), "FACTS_CLIMATE_FORCING.csv.head")
	outfile = os.path.join(os.path.dirname(__file__), "FACTS_CLIMATE_FORCING.csv")
	shutil.copyfile(headfile,outfile)
	WriteToCSV(outfile, np.samps, mode="a")

	# Save the preprocessed data to a pickle
	output = {"SAT": SAT, "Time": Time, "NumTensemble": NumTensemble, "Tlen": Tlen, "scenario": scenario, \
		    "baseyear": baseyear, "infile": fname, "facts_data_file": outfile,  "nsamps": nsamps }
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_preprocess.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the preprocess stage for the emulandice module.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--scenario', help="SSP Emissions scenario", required=True)
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
	parser.add_argument('--baseyear', help="Year to which projections should be referenced", type=int, default=2005)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_preprocess(args.pipeline_id, args.scenario, args.climate_data_file, args.baseyear )

	# Done
	sys.exit()
