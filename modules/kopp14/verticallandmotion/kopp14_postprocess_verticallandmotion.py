import numpy as np
import pickle
import os
import sys
import time
import argparse
import re
from scipy.stats import norm
from read_locationfile import ReadLocationFile
from netCDF4 import Dataset

''' kopp14_postprocess_verticallandmotion.py

This runs the post-processing stage for the vertical land motion component of the Kopp 14
workflow. Projections generated from this stage are site-specific.

Parameters:
nsamps = Number of samples to draw
rng_seed = Seed value for the random number generator
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

'''

def angd(lat0, lon0, qlat, qlon):
	
	# Convert the input from degrees to radians
	(lat0, lon0) = np.radians((lat0, lon0))
	(qlat, qlon) = np.radians((qlat, qlon))
	
	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(qlat)*np.sin(qlon-lon0))**2 + \
	(np.cos(lat0)*np.sin(qlat) - np.sin(lat0)*np.cos(qlat) * np.cos(qlon-lon0))**2),\
	(np.sin(lat0)*np.sin(qlat) + np.cos(lat0)*np.cos(qlat)*np.cos(qlon-lon0)))
	
	# Convert the results from radians to degrees and return
	return(np.degrees(temp))

def NearestPoint(qlat, qlon, lats, lons, tol = None):
	
	# Get the distance between the query point and all the possible points
	dist = angd(lats, lons, qlat, qlon)
	
	# Which is the closest point
	nearest_idx = np.argmin(dist)

	# Is the point within the tolerance?
	if isinstance(tol, (int, float)):
		if dist[nearest_idx] > tol:
			return(None)
	
	return(nearest_idx)


def NearestPoints(qlats, qlons, lats, lons, tol):
	
	if len(qlats) != len(qlons):
		raise Exception("Query lats ({}) and lons ({}) differ in length".format(len(qlats), len(qlons)))
	
	idx = map(lambda qlat,qlon: NearestPoint(qlat, qlon, lats, lons, tol), qlats, qlons)
	
	return(list(idx))
	

def kopp14_postprocess_verticallandmotion(nsamps, rng_seed, locationfilename, pipeline_id):

	# Read in the data from the preprocessing stage
	datafile = "{}_data.pkl".format(pipeline_id)
	try:
		f = open(datafile, 'rb')
	except:
		print("Cannot open datafile\n")
	
	# Extract the data from the file
	my_data = pickle.load(f)
	
	# Extract the relevant data
	names = my_data['names']
	ids = my_data['ids']
	lats = my_data['lats']
	lons = my_data['lons']
	rates = my_data['rates']
	sds = my_data['sds']
	baseyear = my_data['baseyear']
	
	# Define the target years
	targyears = np.arange(2000, 2101, 10)
	
	# Load site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
	# Find the nearest points for the query lats/lons
	site_ids_map = NearestPoints(site_lats, site_lons, lats, lons, tol=None)
	
	# Evenly sample an inverse normal distribution
	np.random.seed(rng_seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = np.random.permutation(norm_inv)
	
	# Output quantiles of interest
	out_q = np.unique(np.append(np.linspace(0,1,101), (0.001, 0.005, 0.01, 0.05, 0.167, 0.50, 0.833, 0.95, 0.99, 0.995, 0.999)))
	nq = len(out_q)
	
	# Initialize variable to hold the samples
	local_sl_q = np.full((nq, len(site_ids), len(targyears)), np.nan)
	
	# Loop over the sites
	nsites = len(site_ids)
	for i in np.arange(0,nsites):
		
		# Skip this site if a match wasn't found
		if site_ids_map[i] is None:
			continue
		
		# This site index
		this_site_ind = site_ids_map[i]
		
		# Loop over the target years
		ntimes = len(targyears)
		for j in np.arange(0,ntimes):
			
			# This target year
			targyear = targyears[j]
		
			# Calculate the samples for this location and time
			GIAproj = rates[this_site_ind] * (targyear - baseyear)
			GIAprojsd = sds[this_site_ind] * (targyear - baseyear)
			these_samps = GIAproj + norm_inv_perm*GIAprojsd
			local_sl_q[:,i,j] = np.quantile(these_samps, out_q)
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", nsites)
	year_dim = rootgrp.createDimension("years", ntimes)
	q_dim = rootgrp.createDimension("quantiles", nq)

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable("localSL_quantiles", "i2", ("quantiles", "nsites", "years"), zlib=True, complevel=4)
	localslq.scale_factor = 0.1
	
	# Assign attributes
	rootgrp.description = "Local SLR contributions from vertical land motion according to Kopp 14 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {}".format(pipeline_id)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	localslq.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	q_var[:] = out_q
	localslq[:,:,:] = local_sl_q

	# Close the netcdf
	rootgrp.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp 14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp 14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the postprocessing stage
	kopp14_postprocess_verticallandmotion(args.nsamps, args.seed, args.locationfile, args.pipeline_id)
	
	# Done
	exit()