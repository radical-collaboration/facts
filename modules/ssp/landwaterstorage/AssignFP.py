import numpy as np
import pickle
from netCDF4 import Dataset
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
	
	# Open the fingerprint file
	try:
		nc_fid = Dataset(fp_filename, 'r')
	except:
		print("Cannot open fingerprint file\n")
		raise
	
	# Read in the fingerprint data
	fp = nc_fid.variables['GROUND'][0,:,:]
	fp_lats = nc_fid.variables['lat'][:]
	fp_lons = nc_fid.variables['lon'][:]
	
	# Interpolate the fingerprint to these locations
	fp_interp = interp.RectBivariateSpline(fp_lats, fp_lons, fp, kx=1, ky=1)
	fp_sites = fp_interp.ev(qlats, qlons+360)/100
	
	return(fp_sites)

if __name__ == '__main__':
	
	infile = "./data/REL_GROUNDWATER_NOMASK.nc"
	
	site_lats = np.array([40.1, 32.6, 40.1, 40.1])
	site_lons = np.array([-72.3, -117.8, -72.3, -72.3])
	
	print(AssignFP(infile, site_lats, site_lons))