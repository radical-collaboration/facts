import numpy as np
import argparse
import pickle
import sys
import os
from netCDF4 import Dataset
import time
from scipy.stats import truncnorm

''' FittedISMIP_project_icesheet.py

Runs the FittedISMIP icesheet projection stage.

Parameters:
nsamps              Number of samples to produce
pyear_start			Projection start year
pyear_end			Projection end year
pyear_step			Stepping from projection year start to end
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def FittedISMIP_project_icesheet(nsamps, pyear_start, pyear_end, pyear_step, cyear_start, cyear_end, baseyear, pipeline_id, rngseed):

	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["years"]
	temp_data = my_data["temp_data"]
	scenario = my_data["scenario"]

	# Load the fit file
	datafilename = "{}_fit.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	groups_dict = my_data["groups_dict"]
	models_dict = my_data["models_dict"]
	betas_dict = my_data["betas_dict"]
	sigmas_dict = my_data["sigmas_dict"]
	trend_mean = my_data["trend_mean"]
	trend_sd = my_data["trend_sd"]

	# Extract the ice sources from the fitted dictionaries
	icesources = betas_dict.keys()

	# Define the target projection years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Find the data years that overlap with the target projection years
	(_, datayr_idx, targyear_idx) = np.intersect1d(years, targyears, return_indices = True)

	# Zero out the temperature data to the base year (Fitted models have 0-forced intercept)
	baseyear_idx = np.flatnonzero(years == baseyear)
	if baseyear_idx.size == 0:
		raise Exception("baseyear is not found in temperature data. baseyear = {}".format(baseyear))
	temp_data = temp_data - temp_data[:,baseyear_idx]

	# Set the seed for the RNG
	rng = np.random.default_rng(rngseed)

	# Initialize the samples dictionary to pass to the post-processing stage
	samps_dict = {}

	# Generate the indices for the temperature samples
	#temp_sample_idx = rng.choice(np.arange(temp_data.shape[0]), nsamps)
	temp_sample_idx = np.arange(nsamps)

	# Generate a list of quantiles for the trend samples
	trend_q = rng.random(nsamps)

	# Loop over the ice sources
	#for icesource in icesources:
	for icesource in ["GIS"]:

		# Calculate the trend contributions over time for this ice sheet component
		ice_trend = truncnorm.ppf(trend_q, a=0.0, b=99999.9, loc=trend_mean[icesource], scale=trend_sd[icesource])[:,np.newaxis] * (targyears - baseyear)[np.newaxis,:]

		# Which model parameters do we need
		betas = betas_dict[icesource]
		sigmas = sigmas_dict[icesource]

		# Generate the indices for the model samples
		model_sample_idx = rng.choice(np.arange(betas.shape[0]), nsamps)

		# Loop over the number of samples we need
		samps = []
		temp_samps = []
		time_samps = []
		const_samps = []

		for tidx, midx in zip(temp_sample_idx, model_sample_idx):

			# Generate a sample
			(this_sample, samp_temp, samp_time, samp_const) = my_model(temp_data[tidx,datayr_idx], \
									betas[midx,:], sigmas[midx], \
									targyears - baseyear, pyear_step,rng)
			samps.append(this_sample)
			#temp_samps.append(samp_temp)
			#time_samps.append(samp_time)
			#const_samps.append(samp_const)

		# Convert the sample array into a numpy array
		samps = np.array(samps)
		#temp_samps = np.array(temp_samps)
		#time_samps = np.array(time_samps)
		#const_samps = np.array(const_samps)

		# Add the trend to the samples
		samps = samps + ice_trend

		# If the user wants to extrapolate projections based on rates, do so here
		if cyear_start or cyear_end:
			for i in np.arange(nsamps):
				samps[i,:] = ExtrapolateRate(samps[i,:], targyears, cyear_start, cyear_end)

		# Add the total samples to the samples dictionary
		samps_dict[icesource] = samps

		# Write the global projections to output netCDF files
		WriteNetCDF(samps, icesource, targyears[targyear_idx], scenario, pipeline_id, baseyear)
		#WriteNetCDF(temp_samps, "{}TEMP".format(icesource), targyears[targyear_idx], scenario, pipeline_id)
		#WriteNetCDF(time_samps, "{}TIME".format(icesource), targyears[targyear_idx], scenario, pipeline_id)
		#WriteNetCDF(const_samps, "{}CONST".format(icesource), targyears[targyear_idx], scenario, pipeline_id)

	# Write the combined AIS projections to netcdf
	#ais_samps = samps_dict["WAIS"] + samps_dict["EAIS"] + samps_dict["PEN"]
	#WriteNetCDF(ais_samps, "AIS", targyears[targyear_idx], scenario, pipeline_id, baseyear)

	# Store the variables in a pickle
	output = {'samps_dict': samps_dict, 'scenario': scenario, 'targyears': targyears[targyear_idx], 'baseyear': baseyear}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(None)



def ExtrapolateRate(sample, targyears, cyear_start, cyear_end):

	# If only one of the constant rate years is provided, imply the other
	if cyear_start and not cyear_end:
		cyear_end = cyear_start + 20
	if cyear_end and not cyear_start:
		cyear_start = cyear_end - 20

	# Find the start and end projection values for the rate calculation
	proj_start = np.interp(cyear_start, targyears, sample)
	proj_end = np.interp(cyear_end, targyears, sample)

	# Calculate the rate
	rate = (proj_end - proj_start) / (cyear_end - cyear_start)

	# Make a new projection
	ext_sample = sample
	ext_sample[targyears >= cyear_end] = proj_end + (rate * (targyears[targyears >= cyear_end] - cyear_end))

	# Return this sample
	return(ext_sample)



def my_model(temp, beta, sigma, dyears, delta_time, rng):

	# If the last temperature value is nan, replace it with a linear extrapolation
	if np.isnan(temp[-1]):
		temp[-1] = temp[-2] + (temp[-2] - temp[-3])

	# Produce a projection for this temperature trajectory
	# NOTE - The fitted rates are per year, so multiply by delta_time (pyear_step)
	# to produce the sample.
	dsle_hat_const = (np.ones(temp.shape) * beta[0])* delta_time
	dsle_hat_temp = (beta[1] * temp + beta[2] * temp**2 + beta[3] * temp**3) * delta_time
	dsle_hat_time = (beta[4] * dyears + beta[5] * dyears**2) * delta_time

	# Sum up the individual changes over time
	sle_hat_const = np.cumsum(dsle_hat_const)
	sle_hat_temp = np.cumsum(dsle_hat_temp)
	sle_hat_time = np.cumsum(dsle_hat_time)
	sle_hat = sle_hat_temp + sle_hat_time + sle_hat_const

	# Apply the error from the fit to this projection
	spread = (sigma * 0.0) / 100.0
	#spread = 0.75
	pct_error = rng.uniform(-spread, spread)
	sle_hat *= 1 + pct_error

	return(sle_hat, sle_hat_temp, sle_hat_time, sle_hat_const)


def WriteNetCDF(icesamps, icetype, data_years, scenario, pipeline_id, baseyear):

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_{1}_globalsl.nc".format(pipeline_id, icetype))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	nyr = len(data_years)
	nsamps = icesamps.shape[0]
	year_dim = rootgrp.createDimension("years", nyr)
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
	rootgrp.description = "Global SLR contribution from {} according to FittedISMIP icesheet workflow".format(icetype)
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}".format(pipeline_id)
	rootgrp.scenario = scenario
	rootgrp.baseyear = baseyear
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = data_years
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = icesamps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()

	return(0)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the FittedISMIP icesheet projection stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 10000)", default=10000, type=int)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--crateyear_start', help="Constant rate calculation for projections starts at this year", default=None, type=int)
	parser.add_argument('--crateyear_end', help="Constant rate calculation for projections ends at this year", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Year to which projections are referenced [default=2005]", default=2005, type=int)
	parser.add_argument('--seed', help="Seed for the random number generator [default = 1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the projection stage with the provided arguments
	FittedISMIP_project_icesheet(args.nsamps, args.pyear_start, args.pyear_end, args.pyear_step, args.crateyear_start, args.crateyear_end, args.baseyear, args.pipeline_id, args.seed)

	sys.exit()