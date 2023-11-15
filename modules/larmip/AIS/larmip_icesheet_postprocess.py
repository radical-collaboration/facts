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


def larmip_postprocess_icesheet(locationfilename, chunksize, pipeline_id):

	# Read in the projections
	projfile = "{}_projections.pkl".format(pipeline_id)
	with open(projfile, 'rb') as f:
		my_data = pickle.load(f)

	# Extract the data from the file
	waissamps = my_data['sl_r2'] + my_data['sl_r3'] + my_data['sl_r4'] + my_data['sl_r5']
	eaissamps = my_data['sl_r1']
	targyears = my_data['targyears']
	scenario = my_data['scenario']
	baseyear = my_data['baseyear']
	f.close()

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

	# Get some dimension data from the loaded data structures
	nsamps = eaissamps.shape[0]
	nyears = len(targyears)
	nsites = len(site_ids)

	# Get the fingerprints for all sites from all ice sheets
	fpdir = os.path.join(os.path.dirname(__file__), "FPRINT")
	waisfp = da.array(AssignFP(os.path.join(fpdir,"fprint_wais.nc"), site_lats, site_lons))
	eaisfp = da.array(AssignFP(os.path.join(fpdir,"fprint_eais.nc"), site_lats, site_lons))

	# Rechunk the fingerprints for memory
	waisfp = waisfp.rechunk(chunksize)
	eaisfp = eaisfp.rechunk(chunksize)

	# Apply the fingerprints to the projections
	waissl = np.multiply.outer(waissamps, waisfp)
	eaissl = np.multiply.outer(eaissamps, eaisfp)

	# Add up the east and west components for AIS total
	aissl = waissl + eaissl

	# Define the missing value for the netCDF files
	nc_missing_value = np.nan #np.iinfo(np.int16).min

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions from icesheets according to LARMIP2 workflow",
			"history": "Created " + time.ctime(time.time()),
			"source": "FACTS: {}".format(pipeline_id),
			"scenario": scenario,
			"baseyear": baseyear}

	wais_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), waissl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	eais_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), eaissl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	ais_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), aissl, {"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	# Write the netcdf output files
	wais_out.to_netcdf("{0}_{1}_localsl.nc".format(pipeline_id, "WAIS"), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
	eais_out.to_netcdf("{0}_{1}_localsl.nc".format(pipeline_id, "EAIS"), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
	ais_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the LARMIP icesheet module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing for the parameters specified from the command line argument
	larmip_postprocess_icesheet(args.locationfile, args.chunksize, args.pipeline_id)

	# Done
	sys.exit()
