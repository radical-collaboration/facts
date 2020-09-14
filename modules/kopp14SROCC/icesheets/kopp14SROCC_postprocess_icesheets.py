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


''' kopp14SROCC_postprocess_icesheets.py

This script runs the ice sheet post-processing task for the Kopp 2014 SROCC workflow. This task
uses the global projections from the 'kopp14SROCC_project_icesheets' script and applies
spatially resolved fingerprints to the ice sheet contribution. The result is a netCDF4
file that contains spatially and temporally resolved samples of ice sheet contributions 
to local sea-level rise

Parameters: 
samptype = Type of samples to use.  One of (arsamps, basamps, hysamps [default])
locationfilename = File that contains points for localization
pipeline_id = Unique identifer for the pipeline running this code

Output: NetCDF file containing local contributions from ice sheets

'''

def kopp14SROCC_postprocess_icesheets(samptype, locationfilename, pipeline_id):
	
	# Read in the fitted parameters from parfile
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open proj file\n")
		sys.exit(1)
	
	# Extract the data from the file
	my_data = pickle.load(f)
	projdata = my_data[samptype]
	targyears = my_data['targyears']
	baseyear = my_data['baseyear']
	f.close() 
	
	# Read in the scenario from the corr file
	projfile = "{}_corr.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open corr file\n")
		sys.exit(1)
	
	# Extract the data from the file
	my_data = pickle.load(f)
	scenario = my_data['scenario']
	f.close() 
	
	# Load the site locations	
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	
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
	writeNetCDF(gissl, pipeline_id, "GIS", targyears, baseyear, scenario, site_lats, site_lons, site_ids)
	writeNetCDF(waissl, pipeline_id, "WAIS", targyears, baseyear, scenario, site_lats, site_lons, site_ids)
	writeNetCDF(eaissl, pipeline_id, "EAIS", targyears, baseyear, scenario, site_lats, site_lons, site_ids)


def writeNetCDF(data, pipeline_id, icesheet_name, targyears, baseyear, scenario, site_lats, site_lons, site_ids):
	
	# Calculate the quantiles
	out_q = np.unique(np.append(np.linspace(0,1,101), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))
	nq = len(out_q)
	local_sl_q = np.nanquantile(data, out_q, axis=0)
	local_sl_q = np.transpose(local_sl_q, (0,2,1))
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{0}_{1}_localsl.nc".format(pipeline_id, icesheet_name)), "w", format="NETCDF4")

	# Define Dimensions
	nsites = len(site_ids)
	nyears = len(targyears)
	nq = len(out_q)
	site_dim = rootgrp.createDimension("nsites", nsites)
	year_dim = rootgrp.createDimension("years", nyears)
	q_dim = rootgrp.createDimension("quantiles", nq)

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable("localSL_quantiles", "i2", ("quantiles", "nsites", "years"), zlib=True, complevel=4)
	#localslq.scale_factor = 0.1

	# Assign attributes
	rootgrp.description = "Local SLR contributions from icesheets according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "SLR Framework: Kopp 2014 workflow, Scenario: {0}, Baseyear: {1}".format(scenario, baseyear)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
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
	
	return(0)
	
if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp14 SROCC SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--samp_type', help="Type of samples to post-process", choices=['hysamps', 'arsamps', 'basamps'], default="hysamps")
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14SROCC_postprocess_icesheets(args.samp_type, args.locationfile, args.pipeline_id)
	
	# Done
	exit()