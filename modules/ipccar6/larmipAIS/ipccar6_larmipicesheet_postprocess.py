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


''' ipccar6_postprocess_larmipicesheet.py

This script runs the ice sheet post-processing task for the ipccar6 larmip workflow. This task
uses the global projections from the 'ipccar6_project_larmipicesheet' script and applies
spatially resolved fingerprints to the ice sheet contribution. The result is a netCDF4
file that contains spatially and temporally resolved samples of ice sheet contributions
to local sea-level rise

Parameters:
locationfile = File that contains points for localization
pipeline_id = Unique identifer for the pipeline running this code

Output: NetCDF file containing local contributions from ice sheets

'''

def ipccar6_postprocess_larmipicesheet(locationfilename, chunksize, pipeline_id):

	# Read in the fitted parameters from parfile
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)

	# Extract the data from the file
	my_data = pickle.load(f)
	waissamps = my_data['wais_samps']
	eaissamps = my_data['eais_samps']
	pensamps = my_data['pen_samps']
	targyears = my_data['years']
	scenario = my_data['scenario']
	baseyear = my_data['baseyear']
	f.close()

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)

	# Merge the PEN samples into WAIS
	waissamps = waissamps + pensamps

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
			"source": "SLR Framework: LARMIP2 workflow",
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
	ais_out.to_netcdf("{0}_{1}_localsl.nc".format(pipeline_id, "AIS"), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the ipccar6 icesheet SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing for the parameters specified from the command line argument
	ipccar6_postprocess_larmipicesheet(args.locationfile, args.chunksize, args.pipeline_id)

	# Done
	exit()
