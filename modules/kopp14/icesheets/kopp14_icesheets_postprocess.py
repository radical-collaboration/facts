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


''' kopp14_postprocess_icesheets.py

This script runs the ice sheet post-processing task for the Kopp 2014 workflow. This task
uses the global projections from the 'kopp14_project_icesheets' script and applies
spatially resolved fingerprints to the ice sheet contribution. The result is a netCDF4
file that contains spatially and temporally resolved samples of ice sheet contributions 
to local sea-level rise

Parameters: 
samptype = Type of samples to use.  One of (arsamps, basamps, hysamps [default])
locationfilename = File that contains points for localization
pipeline_id = Unique identifer for the pipeline running this code

Output: NetCDF file containing local contributions from ice sheets

'''

def kopp14_postprocess_icesheets(samptype, locationfilename, pipeline_id):
	
	# Read in the fitted parameters from parfile
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
	# Extract the data from the file
	my_data = pickle.load(f)
	projdata = my_data[samptype]
	targyears = my_data['targyears']
	f.close()
	
	# Get the fingerprints for all sites from all ice sheets
	fpdir = os.path.join(os.path.dirname(__file__), "FPRINT")
	gisfp = AssignFP(os.path.join(fpdir,"fprint_gis.nc"), site_lats, site_lons)
	waisfp = AssignFP(os.path.join(fpdir,"fprint_wais.nc"), site_lats, site_lons)
	eaisfp = AssignFP(os.path.join(fpdir,"fprint_eais.nc"), site_lats, site_lons)
	
	# Multiply the fingerprints and the projections
	gissl = np.multiply.outer(projdata[:,:,0], gisfp)
	waissl = np.multiply.outer(projdata[:,:,1], waisfp)
	eaissl = np.multiply.outer(projdata[:,:,2], eaisfp)
	
	# Write to netcdf
	writeNetCDF(gissl, pipeline_id, "GIS", targyears, site_lats, site_lons, site_ids)
	writeNetCDF(waissl, pipeline_id, "WAIS", targyears, site_lats, site_lons, site_ids)
	writeNetCDF(eaissl, pipeline_id, "EAIS", targyears, site_lats, site_lons, site_ids)


def writeNetCDF(data, pipeline_id, icesheet_name, targyears, site_lats, site_lons, site_ids):
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{0}_{1}_localsl.nc".format(pipeline_id, icesheet_name)), "w", format="NETCDF4")

	# Define Dimensions
	nsites = len(site_ids)
	ntimes = len(targyears)
	nsamps = data.shape[0]
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
	rootgrp.description = "Local SLR contributions from icesheets according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "SLR Framework: Kopp 2014 workflow"
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	samps.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:]   = site_lats
	lon_var[:]   = site_lons
	loc_var[:]   = site_ids
	year_var[:]  = targyears
	samp_var[:]  = np.arange(nsamps)
	samps[:,:,:] = data

	# Close the netcdf
	rootgrp.close()
	
	return(0)
	
if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--samp_type', help="Type of samples to post-process", choices=['hysamps', 'arsamps', 'basamps'], default="hysamps")
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14_postprocess_icesheets(args.samp_type, args.locationfile, args.pipeline_id)
	
	# Done
	exit()