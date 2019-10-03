import numpy as np
import pickle
import os
import sys
import time
import argparse
import re
from scipy.stats import norm
from netCDF4 import Dataset

''' kopp14_postprocess_verticallandmotion.py

This runs the post-processing stage for the vertical land motion component of the Kopp14
workflow. Projections generated from this stage are site-specific.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
site_ids = ID numbers of the sites of interest
pipeline_id = Unique identifier for the pipeline running this code

'''

def kopp14_project_verticallandmotion(nsamps, rng_seed, site_ids, pipeline_id):

	# Read in the data from the preprocessing stage
	datafile = "{}_data.pkl".format(pipeline_id)
	try:
		f = open(datafile, 'rb')
	except:
		print("Cannot open datafile\n")
	
	# Extract the data from the file
	my_data = pickle.load(f)
	
	# Extract the relevant data
	targyears = my_data['targyears']
	targregions = my_data['targregions']
	targregionnames = my_data['targregionnames']
	rateprojs = my_data['rateprojs']
	rateprojssd = my_data['rateprojssd']
	baseyear = my_data['baseyear']
	
	# Make sure all the requested IDs are available
	missing_ids = np.setdiff1d(site_ids, targregions)
	if(len(missing_ids) != 0):
		missing_ids_string = ",".join(str(this) for this in missing_ids)
		raise Exception("The following IDs are not available: {}".format(missing_ids_string))
	
	# Map the requested site IDs to target regions
	site_ids_map = np.flatnonzero(np.isin(targregions, site_ids))
	site_ids = targregions[site_ids_map]
	
	# Evenly sample an inverse normal distribution
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	
	# Initialize variable to hold the samples
	samps = np.empty((nsamps, len(targyears), len(site_ids)))
	
	# Loop over the sites
	nsites = len(site_ids)
	for i in np.arange(0,nsites):
		
		# This site index
		this_site_ind = site_ids_map[i]
		
		# Reset the RNG seed
		np.random.seed(rng_seed)
		
		# Create a permutation of the inverse normal distribution
		norm_inv_perm = np.random.permutation(norm_inv)
		
		# Loop over the target years
		ntimes = len(targyears)
		for j in np.arange(0,ntimes):
			
			# This target year
			targyear = targyears[j]
		
			# Calculate the samples for this location and time
			GIAproj = rateprojs[this_site_ind] * (targyear - baseyear)
			GIAprojsd = rateprojssd[this_site_ind] * (targyear - baseyear)
			samps[:,j,i] = GIAproj + norm_inv_perm*GIAprojsd
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", nsites)
	year_dim = rootgrp.createDimension("years", ntimes)
	samp_dim = rootgrp.createDimension("samples", nsamps)

	# Populate dimension variables
	#lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	#lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	localsl = rootgrp.createVariable("localSL", "f4", ("samples", "years", "nsites"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from vertical land motion according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {}".format(pipeline_id)
	#lat_var.units = "Degrees North"
	#lon_var.units = "Degrees West"
	id_var.units = "[-]"
	year_var.units = "[-]"
	samp_var.units = "[-]"
	localsl.units = "mm"

	# Put the data into the netcdf variables
	#lat_var[:] = site_lats
	#lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nsamps)
	localsl[:,:,:] = samps

	# Close the netcdf
	rootgrp.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the Kopp14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--site_ids', help="Site ID numbers (from PSMSL database) to make projections for")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Convert the string of site_ids to a list
	site_ids = [int(x) for x in re.split(",\s*", str(args.site_ids))]
	
	# Run the postprocessing stage
	kopp14_project_verticallandmotion(args.nsamps, args.seed, site_ids, args.pipeline_id)
	
	# Done
	exit()