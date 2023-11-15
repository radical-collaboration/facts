import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from netCDF4 import Dataset
from read_locationfile import ReadLocationFile
from AssignFP import AssignFP

import xarray as xr
import dask.array as da


''' ipccar6_postprocess_gmipemuglaciers.py

This script runs the glacier post-processing task for the GMIP2 emulated workflow. This task
uses the global projections from the 'ipccar6_project_gmipemuglaciers' script and applies
spatially resolved fingerprints to the glacier and ice cap contributions. The result is a
netCDF4 file that contains spatially and temporally resolved samples of GIC contributions
to local sea-level rise

Parameters:
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing local contributions from GIC

'''

def ipccar6_postprocess_gmipemuglaciers(locationfilename, chunksize, pipeline_id):

	# Read in the global projections
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)

	# Extract the projection data from the file
	my_data = pickle.load(f)
	gicsamps = my_data["gic_samps"]
	targyears = my_data["years"]
	scenario = my_data["scenario"]
	baseyear = my_data["baseyear"]
	f.close()

	# Load the fingerprint metadata
	fpfile = os.path.join(os.path.dirname(__file__), "fingerprint_region_map.csv")
	fpmap_data = np.genfromtxt(fpfile, dtype=None, names=True, delimiter=',', encoding=None)

	# Extract the data
	fpmapperids = fpmap_data['IceID']
	fpmaps = fpmap_data['FPID']

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

	# Initialize variable to hold the localized projections
	#(nregions, nsamps, nyears) = gicsamps.shape
	gicsamps = np.transpose(gicsamps, (1,0,2))
	(nsamps, nregions, nyears) = gicsamps.shape
	nsites = len(site_ids)
	#local_sl = np.full((nsites, nsamps, nyears), 0.0)
	local_sl = da.zeros((nsamps, nyears, nsites), chunks=(-1,-1,chunksize))

	# Loop through the GIC regions
	for i in np.arange(nregions):

		# Get the fingerprint file name for this region
		fp_idx = np.flatnonzero(fpmapperids == i+1)
		thisRegion = fpmaps[fp_idx][0]

		# Get the fingerprints for these sites from this region
		regionfile = os.path.join(os.path.dirname(__file__), "FPRINT", "fprint_{0}.nc".format(thisRegion))
		regionfp = da.from_array(AssignFP(regionfile, site_lats, site_lons), chunks=chunksize)

		# Multiply the fingerprints and the projections and add them to the running total
		# over the regions
		local_sl += np.multiply.outer(gicsamps[:,i,:], regionfp)

	# Define the missing value for the netCDF files
	nc_missing_value = np.nan #np.iinfo(np.int16).min

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions from glaciers according to GMIP2 emulated workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: AR5 workflow",
			"scenario": scenario,
			"baseyear": baseyear}

	glac_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), local_sl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	glac_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})




	return(None)



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the GMIP2 emulated projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=20]", type=int, default=20)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing for the parameters specified from the command line argument
	ipccar6_postprocess_gmipemuglaciers(args.locationfile, args.chunksize, args.pipeline_id)

	# Done
	exit()
