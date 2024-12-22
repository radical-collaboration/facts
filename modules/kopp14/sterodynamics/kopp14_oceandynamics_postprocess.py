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

''' kopp14_postprocess_oceandynamics.py

This runs the post-processing stage for the ocean dynamics component of the Kopp14
workflow. Projections generated from this stage are site-specific and include both the
ocean dynamics and thermal expansion contributions.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
pipeline_id = Unique identifier for the pipeline running this code

Note that the value of 'nsamps' and 'rng_seed' are shared between the projection stage
and the post-processing stage when run within FACTS.

'''

def kopp14_postprocess_oceandynamics(nsamps, rng_seed, pipeline_id):

	# Read in the configuration -----------------------------------
	infile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)

	# Extract the relevant data
	targyears = my_data['targyears']
	rcp_scenario = my_data['rcp_scenario']
	GCMprobscale = my_data['GCMprobscale']
	maxDOF = my_data['maxDOF']

	# Read in the ZOS data file ------------------------------------
	infile = "{}_ZOS.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)

	focus_site_ids = my_data['focus_site_ids']
	focus_site_lats = my_data['focus_site_lats']
	focus_site_lons = my_data['focus_site_lons']

	# Read in the TE fit data file --------------------------------
	infile = "{}_thermalexp_fit.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)

	ThermExpMean = my_data['ThermExpMean']
	ThermExpStd = my_data['ThermExpStd']
	ThermExpYears = my_data['ThermExpYears']

	# Read in the OD fit data file --------------------------------
	infile = "{}_oceandynamics_fit.pkl".format(pipeline_id)
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
		raise

	# Extract the data from the file
	my_data = pickle.load(f)

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

	tesamps = my_data['thermsamps']

	# Evenly sample an inverse normal distribution
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)

	# Initialize variable to hold the samples
	local_sl = np.empty((nsamps, len(targyears), len(focus_site_ids)))

	# Determine the thermal expansion scale coefficient
	ThermExpScale = norm.ppf(0.95)/norm.ppf(GCMprobscale)

	# Loop over the sites
	nsites = len(focus_site_ids)
	for i in np.arange(0,nsites):

		# This site index
		this_site_ind = i

		# Reset the RNG seed
		rng = np.random.default_rng(rng_seed)

		# Create a permutation of the inverse normal distribution
		norm_inv_perm = rng.permutation(norm_inv)

		# Loop over the target years
		ntimes = len(targyears)
		for j in np.arange(0,ntimes):

			# This target year
			targyear = targyears[j]

			# Find which index targyear is located in the OD and TE years
			od_year_ind = np.argmin(np.abs(OceanDynYears - targyear))
			te_year_ind = np.argmin(np.abs(ThermExpYears - targyear))

			# Calculate the conditional mean and sd
			condmean = OceanDynMean[od_year_ind,this_site_ind] +  OceanDynStd[od_year_ind,this_site_ind] * OceanDynTECorr[od_year_ind,this_site_ind] * (tesamps[:,j]-ThermExpMean[te_year_ind])/ThermExpStd[te_year_ind]
			condstd = (ThermExpScale * OceanDynStd[od_year_ind,this_site_ind]) * np.sqrt(1 - OceanDynTECorr[od_year_ind,this_site_ind]**2)

			# Make the projection
			local_sl[:,j,i] = t.ppf(norm.cdf(norm_inv_perm),OceanDynDOF[od_year_ind,this_site_ind]) * condstd + condmean


	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", ntimes)
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim  = rootgrp.createDimension("locations", nsites)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var  = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var  = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var  = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from ocean dynamics according to K14 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}. ".format(pipeline_id, rcp_scenario)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	samps.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:]  = focus_site_lats
	lon_var[:]  = focus_site_lons
	loc_var[:]  = focus_site_ids
	year_var[:]  = targyears
	samp_var[:]  = np.arange(nsamps)
	samps[:,:,:] = local_sl

	# Close the netcdf
	rootgrp.close()

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp14 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing stage
	kopp14_postprocess_oceandynamics(args.nsamps, args.seed, args.pipeline_id)

	# Done
	exit()
