import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from netCDF4 import Dataset
from read_bkgdrate import read_bkgdrate
from AssignFP import AssignFP


''' kopp14SROCC_postprocess_icesheets.py

This script runs the ice sheet post-processing task for the Kopp 2014 SROCC workflow. This task
uses the global projections from the 'kopp14SROCC_project_icesheets' script and applies
spatially resolved fingerprints to the ice sheet contribution. The result is a netCDF4
file that contains spatially and temporally resolved samples of ice sheet contributions 
to local sea-level rise

Parameters: 
samptype = Type of samples to use.  One of (arsamps, basamps, hysamps [default])
site_ids = Location ids to fingerprint the icesheet contributions (from PSMSL)
pipeline_id = Unique identifer for the pipeline running this code

Output: NetCDF file containing local contributions from ice sheets

'''

def kopp14SROCC_postprocess_icesheets(samptype, focus_site_ids, pipeline_id):
	
	# Read in the fitted parameters from parfile
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)
	
	# Load the site locations	
	ratefile = os.path.join(os.path.dirname(__file__), "bkgdrate.tsv")
	(_, site_ids, site_lats, site_lons) = read_bkgdrate(ratefile, True)
	
	# Test to make sure the list of sites are valid
	if np.any([x >= 0 for x in focus_site_ids]):
		_, _, site_inds = np.intersect1d(focus_site_ids, site_ids, return_indices=True)
		site_ids = site_ids[site_inds]
		site_lats = site_lats[site_inds]
		site_lons = site_lons[site_inds]
	
	# Extract the data from the file
	my_data = pickle.load(f)
	projdata = my_data[samptype]
	targyears = my_data['targyears']
	f.close()
	
	# Get the fingerprints for all sites from all ice sheets
	fpdir = os.path.join(os.path.dirname(__file__), "FPRINT")
	gisfp = AssignFP(os.path.join(fpdir,"fprint_gis.mn"), site_lats, site_lons)
	waisfp = AssignFP(os.path.join(fpdir,"fprint_wais.mn"), site_lats, site_lons)
	eaisfp = AssignFP(os.path.join(fpdir,"fprint_eais.mn"), site_lats, site_lons)
	
	# Multiply the fingerprints and the projections
	gissl = np.multiply.outer(projdata[:,:,0], gisfp)
	waissl = np.multiply.outer(projdata[:,:,1], waisfp)
	eaissl = np.multiply.outer(projdata[:,:,2], eaisfp)
	
	# Write to netcdf
	writeNetCDF(gissl, pipeline_id, "GIS", targyears, site_lats, site_lons, site_ids)
	writeNetCDF(waissl, pipeline_id, "WAIS", targyears, site_lats, site_lons, site_ids)
	writeNetCDF(eaissl, pipeline_id, "EAIS", targyears, site_lats, site_lons, site_ids)


def writeNetCDF(data, pipeline_id, icesheet_name, targyears, site_lats, site_lons, site_ids):
	
	# Calculate the quantiles
	out_q = np.unique(np.append(np.linspace(0,1,101), (0.001, 0.005, 0.01, 0.05, 0.167, 0.5, 0.833, 0.95, 0.99, 0.995, 0.999)))
	nq = len(out_q)
	local_sl_q = np.nanquantile(data, out_q, axis=0)
	local_sl_q = np.transpose(local_sl_q, (0,2,1))
	
	# Calculate the mean and sd of the samples
	local_sl_mean = np.nanmean(data, axis=0).T
	local_sl_sd = np.nanstd(data, axis=0).T	
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{0}_{1}_localsl.nc".format(pipeline_id, icesheet_name)), "w", format="NETCDF4")

	# Define Dimensions
	nsites = local_sl_mean.shape[0]
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
	localslq = rootgrp.createVariable("localSL_quantiles", "f4", ("quantiles", "nsites", "years"), zlib=True, least_significant_digit=2)
	localslmean = rootgrp.createVariable("localSL_mean", "f4", ("nsites", "years"), zlib=True, least_significant_digit=2)
	localslsd = rootgrp.createVariable("localSL_std", "f4", ("nsites", "years"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from icesheets according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "SLR Framework: Kopp 2014 workflow"
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	localslq.units = "mm"
	localslmean.units = "mm"
	localslsd.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	q_var[:] = out_q
	localslq[:,:,:] = local_sl_q
	localslmean[:,:] = local_sl_mean
	localslsd[:,:] = local_sl_sd

	# Close the netcdf
	rootgrp.close()
	
	return(0)
	
if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the Kopp14 SROCC SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--samp_type', help="Type of samples to post-process", choices=['hysamps', 'arsamps', 'basamps'], default="hysamps")
	parser.add_argument('--site_ids', help="Site ID numbers (from PSMSL database) to make projections for")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Convert the string of site_ids to a list
	site_ids = [int(x) for x in re.split(",\s*", str(args.site_ids))]
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14SROCC_postprocess_icesheets(args.samp_type, site_ids, args.pipeline_id)
	
	# Done
	exit()