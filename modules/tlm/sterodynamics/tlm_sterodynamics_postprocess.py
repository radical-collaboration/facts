import numpy as np
import pickle
import os
import sys
import time
import argparse
import re
from scipy.stats import norm
from scipy.stats import t
from netCDF4 import Dataset

import xarray as xr
import dask.array as da


''' tlm_postprocess_oceandynamics.py

This runs the post-processing stage for the ocean dynamics component of the TLM
workflow. Projections generated from this stage are site-specific and include both the
ocean dynamics and thermal expansion contributions.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
pipeline_id = Unique identifier for the pipeline running this code

Note that the value of 'nsamps' and 'rng_seed' are shared between the projection stage
and the post-processing stage when run within FACTS.

'''

def tlm_postprocess_oceandynamics(nsamps, rng_seed, chunksize, keep_temp, pipeline_id):

	# Read in the configuration -----------------------------------
	infile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()

	# Extract the relevant data
	targyears = my_data['targyears']
	scenario = my_data['scenario']
	baseyear = my_data['baseyear']
	GCMprobscale = my_data['GCMprobscale']
	no_correlation = my_data['no_correlation']

	# Read in the ZOS data file ------------------------------------
	infile = "{}_ZOS.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()

	site_ids = my_data['focus_site_ids']
	site_lats = my_data['focus_site_lats']
	site_lons = my_data['focus_site_lons']
	zos_modellist = my_data['zos_modellist']
	zos_scenariolist = my_data['zos_scenariolist']

	# Read in the TE fit data file --------------------------------
	infile = "{}_thermalexp_fit.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()

	ThermExpMean = my_data['ThermExpMean']
	ThermExpStd = my_data['ThermExpStd']
	ThermExpYears = my_data['ThermExpYears']
	ThermExpN = my_data['ThermExpN']

	# Read in the OD fit data file --------------------------------
	infile = "{}_oceandynamics_fit.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()

	OceanDynYears = my_data['OceanDynYears']
	OceanDynMean = my_data['OceanDynMean']
	OceanDynStd = my_data['OceanDynStd']
	OceanDynN = my_data['OceanDynN']
	OceanDynTECorr = my_data['OceanDynTECorr']
	OceanDynDOF = my_data['OceanDynDOF']

	# Read in the TE projections data file ------------------------
	infile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()

	tesamps = my_data['thermsamps']


	# Evenly sample quantile space and permutate
	quantile_samps = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	rng = np.random.default_rng(rng_seed)
	quantile_samps = rng.permutation(quantile_samps)

	# Get the number of years and locations
	nyears = len(targyears)
	nsites = len(site_ids)

	# Determine the thermal expansion scale coefficient
	ThermExpScale = norm.ppf(0.95)/norm.ppf(GCMprobscale)

	# Find which index targyear is located in the OD and TE years
	od_year_ind = np.array([np.argmin(np.abs(OceanDynYears - x)) for x in targyears])[:,np.newaxis]
	te_year_ind = np.array([np.argmin(np.abs(ThermExpYears - x)) for x in targyears])

	# Define the missing value for the netCDF files
	nc_missing_value = np.nan #np.iinfo(np.int16).min

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions from thermal expansion and dynamic sea-level according to Kopp 2014 CMIP6/TLM workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: Kopp 2014 CMIP6/TLM workflow",
			"scenario": scenario,
			"baseyear": baseyear}

	# Loop over the chunks of samples
	for i in np.arange(0,nsites,chunksize):

		# This chunk's temporary file name
		temp_filename = "{0}_tempsamps_{1:05d}.nc".format(pipeline_id,int(i/chunksize))

		# If this temporary file exists, move on to the next chunk
		if os.path.isfile(temp_filename):
			print("{} found, skipping to next chunk".format(temp_filename))
			continue

		# Get the sample indices
		top_idx = np.amin([i+chunksize, nsites])
		site_idx = np.arange(i,top_idx)[np.newaxis,:]

		# Calculate the conditional mean and std dev if correlation is needed
		if no_correlation:
			condmean = OceanDynMean[od_year_ind,site_idx]
			condstd = (ThermExpScale * OceanDynStd[od_year_ind,site_idx])
		else:
			condmean = OceanDynMean[od_year_ind,site_idx] +  OceanDynStd[od_year_ind,site_idx] * OceanDynTECorr[od_year_ind,site_idx] * ((tesamps-np.nanmean(tesamps, axis=0))/np.nanstd(tesamps, axis=0))[:,:,np.newaxis]
			condstd = (ThermExpScale * OceanDynStd[od_year_ind,site_idx]) * np.sqrt(1 - OceanDynTECorr[od_year_ind,site_idx]**2)


		# Generate the samples
		samps = t.ppf(quantile_samps[:,np.newaxis,np.newaxis],OceanDynDOF[np.newaxis,od_year_ind,site_idx]) * condstd + condmean
		samps += tesamps[:,:,np.newaxis]

		# Generate the output xarray
		local_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), samps, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats[site_idx[0,:]]),
							"lon": (("locations"), site_lons[site_idx[0,:]])},
							coords={"years": targyears, "locations": site_ids[site_idx[0,:]], "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

		# Write these samples to a temporary netcdf file
		local_out.to_netcdf(temp_filename, encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Open the temporary data sets
	#combined = xr.open_mfdataset("{0}_tempsamps_*.nc".format(pipeline_id), concat_dim="locations", chunks={"locations":chunksize})
	combined = xr.open_mfdataset("{0}_tempsamps_*.nc".format(pipeline_id), chunks={"locations":chunksize})

	# Write the combined data out to the final netcdf file
	combined.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Remove the temporary files
	if keep_temp == 0:
		for i in np.arange(0,nsites,chunksize):
			os.remove("{0}_tempsamps_{1:05d}.nc".format(pipeline_id,int(i/chunksize)))

	# Produce the intermediate data output netCDF file =========================

	# Produce the included model string
	model_string_pieces = ["{0}-{1}".format(zos_modellist[x], zos_scenariolist[x]) for x in np.arange(len(zos_modellist))]
	model_string = "Models and scenarios included: " + ", ".join(model_string_pieces)

	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_OceanDyn.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", nsites)
	odyear_dim = rootgrp.createDimension("OceanDynYears", len(OceanDynYears))

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	odyear_var = rootgrp.createVariable("OceanDynYears", "i4", ("OceanDynYears",))

	# Create a data variable
	oceandynmean = rootgrp.createVariable("OceanDynMean", "f4", ("nsites", "OceanDynYears"))
	oceandynstd = rootgrp.createVariable("OceanDynStd", "f4", ("nsites", "OceanDynYears"))
	oceandynn = rootgrp.createVariable("OceanDynN", "i4", ("nsites", "OceanDynYears"))
	oceandyntecorr = rootgrp.createVariable("OceanDynTECorr", "f4", ("nsites", "OceanDynYears"))
	oceandyndof = rootgrp.createVariable("OceanDynDOF", "i4", ("nsites", "OceanDynYears"))

	# Assign attributes
	rootgrp.description = "Ocean Dynamics intermediate data for the TLM workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}. ".format(pipeline_id, scenario) + model_string
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	odyear_var[:] = OceanDynYears
	oceandynmean[:,:] = OceanDynMean.T
	oceandynstd[:,:] = OceanDynStd.T
	oceandynn[:,:] = OceanDynN.T
	oceandyntecorr[:,:] = OceanDynTECorr.T
	oceandyndof[:,:] = OceanDynDOF.T

	# Close the netcdf
	rootgrp.close()


	return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the TLM ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--keep_temp', help="Keep the temporary files? 1 = keep, 0 = remove [default=0]", type=int, default=0)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing stage
	tlm_postprocess_oceandynamics(args.nsamps, args.seed, args.chunksize, args.keep_temp, args.pipeline_id)

	# Done
	exit()
