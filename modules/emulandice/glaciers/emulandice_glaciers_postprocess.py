import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from read_locationfile import ReadLocationFile
from AssignFP import AssignFP

import xarray as xr
import dask.array as da

''' emulandice_postprocess_glaciers.py

'''

def emulandice_postprocess_glaciers(locationfilename, chunksize, pipeline_id):

	# Read in the projection data
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)

	# Extract the data from the file
	my_data = pickle.load(f)
	gicsamps = my_data['gic_samps']
	targyears = my_data['targyears']
	scenario = my_data['scenario']
	baseyear = my_data['baseyear']
	preprocess_infile = my_data['preprocess_infile']
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
	gicsamps = np.transpose(gicsamps, (1,0,2))
	(nsamps, nregions, nyears) = gicsamps.shape
	nsites = len(site_ids)
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
	ncvar_attributes = {"description": "Local SLR contributions from glaciers according to emulandice glaciers workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: emulandice workflow",
			"scenario": scenario,
			"baseyear": baseyear,
			"preprocess_infile": preprocess_infile}

	glac_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), local_sl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	glac_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the emulandice glaciers SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing for the parameters specified from the command line argument
	emulandice_postprocess_glaciers(args.locationfile, args.chunksize, args.pipeline_id)

	# Done
	exit()
