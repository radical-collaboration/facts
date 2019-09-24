import numpy as np
import sys
import os
import pickle
import time
import argparse
from netCDF4 import Dataset
from read_annual import read_annual
from AssignFP import AssignFP


''' kopp14_postprocess_glaciers.py

This script runs the glacier post-processing task for the Kopp 2014 workflow. This task
uses the global projections from the 'kopp14_project_glaciers' script and applies
spatially resolved fingerprints to the glacier and ice cap contributions. The result is a 
netCDF4 file that contains spatially and temporally resolved samples of GIC contributions 
to local sea-level rise

Parameters: 
pipeline_id = Unique identifier for the pipeline running this code

Output: NetCDF file containing local contributions from GIC

'''

def kopp14_postprocess_glaciers(pipeline_id):
	
	# Read in the global projections
	projfile = "{}_projections.pkl".format(pipeline_id)
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)
	
	# Extract the projection data from the file
	my_data = pickle.load(f)
	gicsamps = my_data["gicsamps"]
	f.close()
	
	# Read in the fingerprint information
	fpfile = "{}_fp.pkl".format(pipeline_id)
	try:
		f = open(fpfile, 'rb')
	except:
		print("Cannot open fpfile\n")
		sys.exit(1)
	
	# Extract the fingerprint information from the file
	my_data = pickle.load(f)
	fpmapperids = my_data["fpmapperids"]
	fpmaps = my_data["fpmaps"]
	f.close()
	
	# Read in the configuration information
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configfile\n")
		sys.exit(1)
	
	# Extract the fingerprint information from the file
	my_data = pickle.load(f)
	targyears = my_data["targyears"]
	rcp_scenario = my_data["rcp_scenario"]
	f.close()
	
	
	# Load the site locations
	rlrdir = os.path.join(os.path.dirname(__file__), "rlr_annual")
	sites = read_annual(rlrdir, True)
	
	# Extract site lats, lons, and ids
	def decomp_site(site):
		return(site.lat, site.lon, site.id)
	vdecomp_site = np.vectorize(decomp_site)
	(site_lats, site_lons, site_ids) = vdecomp_site(sites)
	
	# FOR SIMPLICITY, LOCALIZE TO ONLY A FEW LOCATIONS
	#sites_include = np.array([12,299,396,188,161,10,405,155,43,269,860,526,235,88,1])
	sites_include = np.array([12,299])
	_, _, site_inds = np.intersect1d(sites_include, site_ids, return_indices=True)
	site_ids = site_ids[site_inds]
	site_lats = site_lats[site_inds]
	site_lons = site_lons[site_inds]
	
	# Initialize variable to hold the localized projections
	(nsamps, nregions, ntimes) = gicsamps.shape
	nsites = len(site_ids)
	local_sl = np.full((nsites, nsamps, nregions, ntimes), np.nan) 
	
	# Loop through the GIC regions
	for i in np.arange(0,nregions):
		
		# Get the fingerprint file name for this region
		thisRegion = fpmaps[i]
		
		# Get the fingerprints for these sites from this region
		regionfile = os.path.join(os.path.dirname(__file__), "FPRINT", "fprint_{0}.mn".format(thisRegion))
		regionfp = AssignFP(regionfile, site_lats, site_lons)
		
		# Multiply the fingerprints and the projections
		local_sl[:,:,i,:] = np.transpose(np.multiply.outer(gicsamps[:,i,:], regionfp), (2,0,1))
		
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "{}_localsl.nc".format(pipeline_id)), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", nsites)
	year_dim = rootgrp.createDimension("years", ntimes)
	samp_dim = rootgrp.createDimension("samples", nsamps)
	reg_dim = rootgrp.createDimension("regions", nregions)

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	localsl = rootgrp.createVariable("localSL", "f4", ("nsites", "samples", "regions", "years"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from glaciers and ice caps according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: Kopp 2014 GIC workflow - {0}".format(rcp_scenario)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	id_var.units = "[-]"
	year_var.units = "[-]"
	samp_var.units = "[-]"
	localsl.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nsamps)
	localsl[:,:,:,:] = local_sl

	# Close the netcdf
	rootgrp.close()
	
if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the GIC component of the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected	
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
		
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the postprocessing for the parameters specified from the command line argument
	kopp14_postprocess_glaciers(args.pipeline_id)
	
	# Done
	exit()