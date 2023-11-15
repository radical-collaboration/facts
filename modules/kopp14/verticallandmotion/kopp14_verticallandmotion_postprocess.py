import numpy as np
import pickle
import os
import sys
import time
import argparse
import re
from scipy.stats import norm
from read_locationfile import ReadLocationFile

import xarray as xr
import dask.array as da

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


def kopp14_postprocess_verticallandmotion(nsamps, rng_seed, baseyear, pyear_start, pyear_end, pyear_step, locationfilename, chunksize, pipeline_id):

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

	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	targyears = np.union1d(targyears, baseyear)


	# Load site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

	# Dimension variables
	nyears = len(targyears)
	nsites = len(site_ids)


	# Find the nearest points for the query lats/lons
	site_ids_map = np.array(NearestPoints(site_lats, site_lons, lats, lons, tol=None))

	# Evenly sample an inverse normal distribution
	np.random.seed(rng_seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = np.random.permutation(norm_inv)

	# Missing value for netcdf file
	nc_missing_value = np.iinfo(np.int16).min

	# Get the rates and sds for the locations of interest
	site_rates = da.array([rates[x] for x in site_ids_map])
	site_sds = da.array([sds[x] for x in site_ids_map])

	# Rechunk the rates and sds
	site_rates = site_rates.rechunk(chunksize)
	site_sds = site_sds.rechunk(chunksize)

	# Generate the projected means and standard deviations
	GIAproj = np.multiply.outer(targyears - baseyear, site_rates)
	GIAprojsd = np.multiply.outer(targyears - baseyear, site_sds)

	# Produce the samples from the means and standard deviations
	local_sl = GIAproj + np.multiply.outer(norm_inv_perm, GIAprojsd)

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions from vertical land motion according to Kopp 2014 workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: Kopp 2014 workflow",
			"scenario": "NA",
			"baseyear": baseyear}

	vlm_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), local_sl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	# Write the netcdf output file
	vlm_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})


	return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp 14 vertical land motion workflow",\
	epilog="Note: This is meant to be run as part of the Kopp 14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--baseyear', help="Base or reference year for projetions [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Make sure the base year and target years are within data limits for this module
	if(args.baseyear < 2000):
		raise Exception("Base year cannot be less than year 2000: baseyear = {}".format(args.baseyear))
	if(args.baseyear > 2300):
		raise Exception("Base year cannot be greater than year 2300: baseyear = {}".format(args.baseyear))
	if(args.pyear_start < 2000):
		raise Exception("Projection year cannot be less than year 2000: pyear_start = {}".format(args.pyear_start))
	if(args.pyear_end > 2300):
		raise Exception("Projection year cannot be greater than year 2300: pyear_end = {}".format(args.pyear_end))

	# Make sure the target year stepping is positive
	if(args.pyear_step < 1):
		raise Exception("Projection year step must be greater than 0: pyear_step = {}".format(args.pyear_step))

	# Run the postprocessing stage
	kopp14_postprocess_verticallandmotion(args.nsamps, args.seed, args.baseyear, args.pyear_start, args.pyear_end, args.pyear_step, args.locationfile, args.chunksize, args.pipeline_id)

	# Done
	exit()
