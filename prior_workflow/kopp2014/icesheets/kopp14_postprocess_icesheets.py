import numpy as np
import sys
import os
import pickle
import time
from read_annual import read_annual
from AssignFP import AssignFP
from netCDF4 import Dataset

''' kopp14_postprocess_icesheets.py

This script runs the ice sheet post-processing task for the Kopp 2014 workflow. This task
uses the global projections from the 'kopp14_project_icesheets' script and applies
spatially resolved fingerprints to the ice sheet contribution. The result is a netCDF4
file that contains spatially and temporally resolved samples of ice sheet contributions 
to local sea-level rise

Parameters: 
projfile = Name of the file containing the fitted distribution data
rlrdir = PSMSL data directory for loading site locations
samptype = Type of samples to use.  One of (arsamps, basamps, hysamps [default])

Output: NetCDF file containing local contributions from ice sheets

'''

def kopp14_postprocess_icesheets(projfile, samptype="hysamps"):
	
	# Read in the fitted parameters from parfile
	try:
		f = open(projfile, 'rb')
	except:
		print("Cannot open projfile\n")
		sys.exit(1)
	
	# Load the site locations
	rlrdir = os.path.join(os.path.dirname(__file__), "data", "rlr_annual")
	sites = read_annual(rlrdir, True)
	
	# Extract site lats, lons, and ids
	def decomp_site(site):
		return(site.lat, site.lon, site.id)
	vdecomp_site = np.vectorize(decomp_site)
	(site_lats, site_lons, site_ids) = vdecomp_site(sites)
	
	# Extract the data from the file
	my_data = pickle.load(f)
	projdata = my_data[samptype]
	f.close()
	
	# Get the fingerprints for all sites from all ice sheets
	fpdir = os.path.join(os.path.dirname(__file__), "data", "FPRINT")
	gisfp = AssignFP(os.path.join(fpdir,"fprint_gis.mn"), site_lats, site_lons)
	waisfp = AssignFP(os.path.join(fpdir,"fprint_wais.mn"), site_lats, site_lons)
	eaisfp = AssignFP(os.path.join(fpdir,"fprint_eais.mn"), site_lats, site_lons)
	
	# Multiply the fingerprints and the projections
	gissl = np.multiply.outer(projdata[:,:,0], gisfp)
	waissl = np.multiply.outer(projdata[:,:,1], waisfp)
	eaissl = np.multiply.outer(projdata[:,:,2], eaisfp)
	
	# Write the localized projections to a netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "data", "local_slr.nc"), "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", len(sites))
	year_dim = rootgrp.createDimension("years", projdata.shape[1])
	samp_dim = rootgrp.createDimension("samples", projdata.shape[0])

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	localgis = rootgrp.createVariable("localGIS", "f4", ("samples", "years", "nsites"), zlib=True, least_significant_digit=2)
	localwais = rootgrp.createVariable("localWAIS", "f4", ("samples", "years", "nsites"), zlib=True, least_significant_digit=2)
	localeais = rootgrp.createVariable("localEAIS", "f4", ("samples", "years", "nsites"), zlib=True, least_significant_digit=2)

	# Assign attributes
	rootgrp.description = "Local SLR contributions from icesheets according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "SLR Framework: Kopp 2014 workflow"
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees West"
	id_var.units = "[-]"
	year_var.units = "[-]"
	samp_var.units = "[-]"
	localgis.units = "mm"
	localwais.units = "mm"
	localeais.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = site_lats
	lon_var[:] = site_lons
	id_var[:] = site_ids
	year_var[:] = 2010 + 10*(np.arange(0,projdata.shape[1]))
	samp_var[:] = np.arange(0,projdata.shape[0])
	localgis[:,:,:] = gissl
	localwais[:,:,:] = waissl
	localeais[:,:,:] = eaissl

	# Close the netcdf
	rootgrp.close()
	
if __name__ == '__main__':
	
	# Run the postprocessing for the parameters specified from the command line argument
	try:
		projfile = sys.argv[1]
	except:
		projfile = os.path.join(os.path.dirname(__file__), "data", "kopp14_icesheets_projections.pkl")
		
	try:
		samptype = sys.argv[2]
	except:
		samptype = "hysamps"
	finally:
		kopp14_postprocess_icesheets(projfile, samptype)
	
	# Done
	exit()