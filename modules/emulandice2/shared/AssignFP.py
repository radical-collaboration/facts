import numpy as np
from ReadFingerprint import ReadFingerprint as readfp
from scipy import interpolate as interp

''' AssignFP.py

Assigns interpolated fingerprint coefficients to sites identified by the vectors
of lats and lons provided.

Parameters:
fp_filename = Fingerprint file passed to ReadFingerprint
qlats = Vector of latitudes of sites of interest [-90, 90]
qlons = Vector of longitudes of sites of interest [-180, 180]

Return:
fp_sites = Vector of fingerprint coefficients for the sites of interest


'''

def AssignFP(fp_filename, qlats, qlons):
	## Read in the fitted parameters from parfile
	# Open the file
	try:
		(fp, fp_lats, fp_lons) = readfp(fp_filename)
	except:
		print("Cannot open fingerprint file\n")

	# Interpolate the fingerprint to these locations
	lat_sort = np.argsort(fp_lats)
	fp_interp = interp.RectBivariateSpline(fp_lats[lat_sort], fp_lons, fp[lat_sort,:], kx=1, ky=1)
	fp_sites = fp_interp.ev(qlats, np.mod(qlons, 360))*1000
	
	return(fp_sites)
