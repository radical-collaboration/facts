import numpy as np
import os
import sys
import fnmatch
import argparse
import re
import pickle
import time
import subprocess
import shutil
from netCDF4 import Dataset
from scipy.stats import truncnorm


# For AIS, there are three regions (WAIS, EAIS, and PEN)
def ExtractProjections(emulandice_file):

	# Initialize
	ice_sources = []
	regions = []
	years = []
	samples = []
	sles = []

	# Open the emulandice_file for reading
	with open(emulandice_file, "r") as f:

		# Skip the header line
		_ = f.readline()

		# Load the rest of the data
		for line in f:
			lp = re.split(",", line)
			ice_sources.append(lp[0])
			regions.append(lp[1])
			years.append(int(lp[2]))
			samples.append(int(lp[3]))
			sles.append(float(lp[7]))

	# Convert to numpy arrays
	ice_sources = np.array(ice_sources)
	regions = np.array(regions)
	years = np.array(years)
	samples = np.array(samples)
	sles = np.array(sles)

	# Extract the target years
	targyears = np.unique(years)
	unique_samples = np.unique(samples)

	# Initialize the return data structure
	wais_data = np.full((len(unique_samples), len(targyears)), np.nan)
	eais_data = np.full((len(unique_samples), len(targyears)), np.nan)
	pen_data = np.full((len(unique_samples), len(targyears)), np.nan)

	# Loop over all the entries
	for i in np.arange(len(sles)):

		# This sample
		this_icesource = ice_sources[i]
		this_region = regions[i]
		this_year = years[i]
		this_sample = samples[i]
		this_sle = sles[i]

		# Skip this entry if it's not for AIS
		if this_icesource != "AIS":
			continue

		# Which index will this entry be filling?
		year_idx = np.flatnonzero(targyears == this_year)
		sample_idx = this_sample - 1

		# Put the data into the return data structure
		if this_region == "WAIS":
			wais_data[sample_idx, year_idx] = this_sle * 10.0  # Convert cm to mm
		elif this_region == "EAIS":
			eais_data[sample_idx, year_idx] = this_sle * 10.0  # Convert cm to mm
		else:
			pen_data[sample_idx, year_idx] = this_sle * 10.0  # Convert cm to mm

	# Done
	return(wais_data, eais_data, pen_data, targyears)


def emulandice_project_AIS(pipeline_id, icesource="AIS"):

	# Load the preprocessed data
	preprocess_file = "{}_preprocess.pkl".format(pipeline_id)
	with open(preprocess_file, 'rb') as f:
		preprocess_data = pickle.load(f)

	preprocess_infile = preprocess_data["infile"]
	baseyear = preprocess_data["baseyear"]
	scenario = preprocess_data["scenario"]
	nsamps = preprocess_data["nsamps"]

	# Load the fit data
	fit_file = "{}_fit.pkl".format(pipeline_id)
	with open(fit_file, 'rb') as f:
		fit_data = pickle.load(f)

	trend_mean = fit_data["trend_mean"]
	trend_sd = fit_data["trend_sd"]

	# Run the module using the FACTS forcing data

	py_working_dir = os.getcwd()
	emulandice_dataset = 'FACTS_CLIMATE_FORCING.csv'
	subprocess.run(["bash", "emulandice_steer.sh", emulandice_dataset, str(nsamps), icesource])

	# Get the output from the emulandice run
	emulandice_file = os.path.join(os.path.dirname(__file__),"results", "projections_FAIR_FACTS.csv")
	wais_samples, eais_samples, pen_samples, targyears = ExtractProjections(emulandice_file)

	# Make sure we get the number of samples we expected
	if nsamps != wais_samples.shape[0]:
		raise Exception("Number of SLC projections does not match number of temperature trajectories: {} != {}".format(wais_samples.shape[0], nsamps))

	# Generate samples for trends correlated among ice sheets
	# Note: Keep seed hard-coded and matched with GrIS module within emulandice module set
	rng = np.random.default_rng(8071)
	trend_q = rng.random(nsamps)

	# Calculate the trend contributions over time for each ice sheet component
	eais_trend = truncnorm.ppf(trend_q, a=0, b=99999, loc=trend_mean["EAIS"], scale=trend_sd["EAIS"])[:,np.newaxis] * (targyears - baseyear)[np.newaxis,:]
	wais_trend = truncnorm.ppf(trend_q, a=0, b=99999, loc=trend_mean["WAIS"], scale=trend_sd["WAIS"])[:,np.newaxis] * (targyears - baseyear)[np.newaxis,:]
	pen_trend = truncnorm.ppf(trend_q, a=0, b=99999, loc=trend_mean["PEN"], scale=trend_sd["PEN"])[:,np.newaxis] * (targyears - baseyear)[np.newaxis,:]

	# Store the samples for AIS components
	eais_samples += eais_trend
	wais_samples += wais_trend
	pen_samples += pen_trend

	# Save the global projections to a pickle
	output = {"waissamps": wais_samples + pen_samples, "eaissamps": eais_samples, "targyears": targyears, "scenario": scenario, \
			"baseyear": baseyear, "preprocess_infile": preprocess_infile}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the global projections to netcdf files
	WriteNetCDF(eais_samples+wais_samples+pen_samples, None, targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(eais_samples, "EAIS", targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(pen_samples, "PEN", targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(wais_samples, "WAIS", targyears, baseyear, scenario, nsamps, pipeline_id)

	# Done
	return(None)


def WriteNetCDF(slr, region, targyears, baseyear, scenario, nsamps, pipeline_id):

	# Write the total global projections to a netcdf file
	if region is None:
		nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
		nc_description = "Global SLR contribution from Antarctica using the emulandice module"
	else:
		nc_filename = os.path.join(os.path.dirname(__file__), "{}_{}_globalsl.nc".format(pipeline_id, region))
		nc_description = "Global SLR contribution from Antarctica ({}) using the emulandice module".format(region)
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = nc_description
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}. ".format(pipeline_id)
	rootgrp.baseyear = baseyear
	rootgrp.scenario = scenario
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = slr[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	return(None)



if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the emulandice AIS SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_project_AIS(args.pipeline_id)

	# Done
	sys.exit()
