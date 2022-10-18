import numpy as np
import sys
from netCDF4 import Dataset

''' ReadFingerprint.py

Provides a function that reads in a fingerprint data file from the netCDF files created
from the original spherical harmonics files.

Parameters:
fname = File name that contains the lat/lon and fingerprint information

Return:
f = The fingerprint coefficient along a lat/lon grid [nlat x nlon]
lat = Vector of latitudes
lon = Vector of longitudes

'''

def ReadFingerprint(fname):

	# Open the fingerprint file
	try:
		nc_fid = Dataset(fname, 'r')
	except:
		print("Cannot open fingerprint file: {0}\n".format(fname))
		raise
	
	# Read in the fingerprint data
	fp = nc_fid.variables['fp'][:,:]
	fp_lats = nc_fid.variables['lat'][:]
	fp_lons = nc_fid.variables['lon'][:]

	return(fp, fp_lats, fp_lons)
	