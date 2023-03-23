import numpy as np
import pickle
import sys
import os
import argparse
import time
import re
from read_locationfile import ReadLocationFile

import xarray as xr
import dask.array as da

''' kopp14_postprocess_landwaterstorage.py

This script runs the thermal expansion postprocessing task for the Kopp 2014 workflow.
This task generates localized contributions to sea-level change due to land water storage.

Parameters:
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing the local sea-level rise projections

Note: The Kopp14 workflow does not apply a fingerprint to the land water storage component.
As such, this script applies an implicit fingerprint coefficient of '1' and produces
localized slr projections for the site ids provided.

'''

def kopp14_postprocess_landwaterstorage(locationfilename, chunksize, pipeline_id):

	# Load the configuration file
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projection file {}\n").format(projfile)

	# Extract the configuration variables
	my_proj = pickle.load(f)
	f.close()

	targyears = my_proj["years"]
	lwssamps = np.transpose(my_proj["lwssamps"])
	baseyear = my_proj["baseyear"]

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

	# Initialize variable to hold the localized projections
	nsamps = lwssamps.shape[0]
	nyears = len(targyears)
	nsites = len(site_ids)

	# Apply the effective fingerprint of "1" to the global projections for each site
	fpsites = da.ones((nsites), chunks=chunksize)
	local_sl = np.multiply.outer(lwssamps, fpsites)

	# Define the missing value for the netCDF files
	nc_missing_value = np.iinfo(np.int16).min

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions from land water storage according to Kopp 2014 workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: Kopp 2014 workflow",
			"scenario": "NA",
			"baseyear": baseyear}

	lws_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), local_sl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	lws_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})


	return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the land water storage postprocessing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the projection process on the files specified from the command line argument
	kopp14_postprocess_landwaterstorage(args.locationfile, args.chunksize, args.pipeline_id)

	exit()
