import numpy as np
import os
import iris
import sys
import argparse
import fnmatch
import re


'''
Apply an inverse-barometer-effect correction to ZOS data.

See below paper for details:
Piecuch, C. G., and Ponte, R. M. (2015), 
Inverted barometer contributions to recent sea level changes along the northeast coast 
of North America, Geophys. Res. Lett., 42, 5918–5925, doi:10.1002/2015GL064580. 

Note: In order to run this code, activate the "esmpy" conda environment. Required for 
access to the "iris" module.

	conda activate esmpy

'''

def guess_bounds(cube, coords):
	"""Guess bounds of a cube, or not."""
	# check for bounds just in case
	for coord in coords:
		if not cube.coord(coord).has_bounds():
			cube.coord(coord).guess_bounds()
	return cube


def calc_area_weights(cube):
	
	# Calculate area weights from the provided grid
	cube = guess_bounds(cube, ['longitude', 'latitude'])
	areaweights = iris.analysis.cartography.area_weights(cube)
	
	# Apply the PSL mask to the area weights
	areaweights = np.ma.masked_where(cube.data.mask, areaweights)
	
	# Normalize the area weights
	norm_area = areaweights / np.sum(areaweights, axis=(1,2))[:,None,None]
	
	return(norm_area[0,...][None,...])
	

def gravity(latitudes):
	
	# WGS (World Geodetic System) 84 Ellipsoidal Gravity Formula
	# https://en.wikipedia.org/wiki/Gravity_of_Earth#Latitude_model
	
	# Define constants
	a = 6378137.0		# Equatorial semi-axis (m)
	b = 6356752.3		# Polar semi-axis (m)
	ge = 9.780327		# Equatorial gravitational acceleration (m s-2)
	gp = 9.8321849378	# Polar gravitational acceleration (m s-2)
	
	# Derive other constants
	e_squared = 1.0 - (b/a)**2	# Spheroid's eccentricity
	k = (b*gp - a*ge) / (a*ge)	# Formula constant
	
	# Convert latitudes from degrees to radians
	lats = latitudes * np.pi/180.0
	
	# Calculate the latitude-dependent gravity
	g_num = 1 + k*np.sin(lats)**2
	g_denom = np.sqrt(1 - e_squared * np.sin(lats)**2)
	g_lat = ge * (g_num/g_denom)
	
	# Return gravity
	return(g_lat)
	

def main(zosdir, psldir, outdir):
	
    # Flag for calculating area weights
    calc_weights_flag = True
	
    # Loop over the models available in the ZOS directory
    #for model in os.listdir(zosdir):
    #for model in (model for model in os.listdir(zosdir) if model not in ['.DS_Store']):
    for model in (model for model in os.listdir(zosdir) if model=='UKESM1-0-LL'):
        # Skip if this isn't really a model
        if re.search(r"^\.", model):
            continue
		
        # Define this ZOS and PSL model directories
        zos_model_dir = os.path.join(zosdir, model)
        psl_model_dir = os.path.join(psldir, model)
		
        # DEBUG: Which model are we on?
        print(model)
		
        # Reset the area weight flag if you want to calculate new area weights for each model
        # Note: Comment this out if you want to assume the same area weights across all models
        calc_weights_flag = True

        # Loop over the files available for this model
        for zos_file in sorted(os.listdir(zos_model_dir)):

			
            # Extract the pieces of the file name
            pieces = zos_file.split("_")
            this_model = pieces[2]
            this_exp = pieces[3]
            this_variant = pieces[4]
            this_grid = pieces[5]
            this_time = pieces[6]
			
			# Is there a matching model available in the PSL data?
            try:
                psl_files = os.listdir(psl_model_dir)
            except:
                print("Cannot open PSL model directory: {}. Moving on...".format(psl_model_dir))
                continue
			
			# Is there a matching PSL file for this model?
            match_psl_file = fnmatch.filter(psl_files, "psl_Amon_{0}_{1}_{2}_*_{3}*.nc".format(this_model, this_exp, this_variant, this_time))
            if not match_psl_file:
                print("Could not find matching PSL file for {}. Moving on...".format(zos_file))
                continue
			
			# A match exists. Load the cubes.
            cube_zos = iris.load_cube(os.path.join(zos_model_dir, zos_file))
            cube_psl = iris.load_cube(os.path.join(psl_model_dir, match_psl_file[0]))
			
			# Find the ZOS mask and apply to the PSL data
            masked_psl = np.ma.masked_where(np.isnan(cube_zos.data), cube_psl.data)
            cube_psl = cube_psl.copy(data=np.float32(masked_psl))
			
			# Calculate the normalized area weights if needed
            if calc_weights_flag:
                norm_area = calc_area_weights(cube_psl)
                calc_weights_flag = False
	
			# Remove the area-weighted mean from the psl data
            weighted_cube = cube_psl * norm_area
            psl_mean_cube = weighted_cube.collapsed(["latitude", "longitude"], iris.analysis.SUM)
            psl_anom_cube = cube_psl - psl_mean_cube
	
			# Calculate the IBE correction
            rho = 1025.0
            g = 9.80665
			#g = gravity(cube_psl.coord("latitude").points)[None,:,None]	# gravity as function of latitude
            ibe = -1 * (psl_anom_cube.data) / (rho*g)
            
			# Apply the correction to the ZOS field
            cube_zos_ibe = cube_zos + ibe
	
			# Make a copy of the original zos cube, update it with the new zos data, and write it to netcdf
            if not os.path.exists(os.path.join(outdir, model)):
                os.mkdir(os.path.join(outdir, model))
            outfile = os.path.join(outdir, model, re.sub(r"\.nc", "_ibe.nc", zos_file))
            #iris.save(cube_zos.copy(data=np.float32(cube_zos_ibe.data)), outfile, zlib=True, complevel=4,fill_value=np.nan)
	
	# Done
    return(None)

if __name__ == "__main__":
	
    '''
	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Apply the Inverse Barometer Effect Correction to ZOS data")
	
	# Add arguments
	parser.add_argument('--zosdir', help="Directory containing the ZOS model-separated data", required=True)
	parser.add_argument('--psldir', help="Directory containing the PSL model-separated data", required=True)
	parser.add_argument('--outdir', help="Output directory", required=True)
	
	# Parse the arguments
	args = parser.parse_args()
	'''
	# Run the code (args.psldir)
    main(<zos_dir>,<psl_dir>,<output_dir>)
	```


