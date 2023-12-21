import numpy as np
import sys
import os
import pickle
import time
import argparse
import re
from read_locationfile import ReadLocationFile
from AssignFP import AssignFP
import yaml

import xarray as xr
import dask.array as da

''' emulandice_postprocess.py

'''

### THIS SCRIPT SHOULD BE MADE MORE GENERIC SO IT
### READS IN A BUNCH OF _global.nc files, figures out which
### FINGERPRINTS GOES WITH WHICH, AND THEN WRITES OUT A
### _localsl.nc FILE. Also needs to appropriately batch
### the ice_source projections together, so eg there is one
### AIS_localsl.nc file produced from scaling the different sectoral
### _globalsl.nc files

'''
PROPOSED APPROACH:

grdfingerprint.yml will be organized in three levels --
top levels is ice sources to which results for localsl will be aggregated
(typically AIS, GrIS, glaciers, but one might imagine wanting to substitute
a version where WAIS, EAIS and AP were broken out).
Second level is region ids, which will be matched with region ids specified
in each netcdf file. Third level is variables, including one describing the 
filename for the fingerprint, eg

- LWS
	- LWS.ALL
 		- fingerprint: fprint_groundwater.nc
- AIS
	- AIS.EAIS
 		- fingerprint: fprint_eais.nc
   	- AIS.WAIS
    		- fingerprint: fprint_wais.nc
      	- AIS.AP
       		- fingerprint: fprint_ap.nc
- GrIS
	- GIS.ALL
 		- fingerprint: fprint_gis.nc
- glaciers
	- GLA.RGI01
 		- fingerprint: fprint_glac1.nc
	- GLA.RGI02
 		- fingerprint: fprint_glac2.nc

and so forth.

emulandice_postprocess takes as an input a list of nc filenames to process
for globalsl files. It loops through the file, reading attributes
that should be added to the file (and/or trying to make sense of filenames)
to assign each file to a ice source and region. Then for each ice source
for which there is at least one file, we go through all those files, 
scaling each one by the appropriate fingerprint file, and produce
the aggregated localsl file.

This is similar to the old emulandice_postprocess_glaciers script,
except we are working with a collection of netcdf files rather than 
a pickle file, and the approach is more generic (such that the same
script can be used for any source region).
'''


def emulandice_postprocess(locationfilename, chunksize, pipeline_id,ncfiles,grdfingerprintfile, scenario, baseyear, targyears):

	matching_ncfiles_dict, grdfingerprint_dict = process_ncfiles_and_grdfingerprintfile(ncfiles, grdfingerprintfile)

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, site_ids, site_lats, site_lons) = ReadLocationFile(locationfile)
	nsites = len(site_ids)

	scaled_ncfiles_dict, targyears = scale_ncfiles_by_fingerprint(matching_ncfiles_dict, grdfingerprint_dict, 
													site_lats, site_lons, chunksize)

	# loop over scaled_ncfiles_dict. For each ice source that has
	# a non-zero array, write out a netcdf file with the localized projections
	for ice_source, scaled_ncfiles in scaled_ncfiles_dict.items():
		if scaled_ncfiles:
			write_localized_projections(scaled_ncfiles, site_ids, site_lats, site_lons, pipeline_id,
							   scenario, baseyear, ice_source, targyears):
	return(None)

'''
process_ncfiles_and_grdfingerprintfile takes in a list of ncfiles and a grdfingerprintfile.
The grdfingerprintfile is a yaml file with a three-level hierarchy. At the
top level is an arbitarily labeled ice source (e.g., AIS, GrIS, glaciers).
At the next level is a set of regions (e.g., AIS.EAIS, AIS.WAIS, glaciers.RGI01).
Each netcdf file may have a 'region' global attribute that matches one of the
defined regions. At the third level is the variable 'fingerprint', which identifies
the name of the fingerprint file to use for that region.

The function should loop through the grdfingerprintfile, and for each region
identify whether one of the netcdf files matches the region. If so, it should
add that to a dictionary that has keys of region names and values that are arrays
of matching netcdf files. The function should return this dictionary and the grdfingerprintfile dicitonary.
'''

def process_ncfiles_and_grdfingerprintfile(ncfiles, grdfingerprintfile):
	# Load the grdfingerprintfile
	with open(grdfingerprintfile, 'r') as file:
		grdfingerprint_dict = yaml.safe_load(file)

	# Initialize the dictionary to store the matching netcdf files
	matching_ncfiles_dict = {}

	# Loop through the grdfingerprintfile
	for ice_source, regions in grdfingerprint_dict.items():
		for region, fingerprint in regions.items():
			# Initialize an empty list for this region
			matching_ncfiles_dict[region] = []

			# Loop through the ncfiles
			for ncfile in ncfiles:
				# Open the netcdf file
				ds = xr.open_dataset(ncfile)

				# Check if the 'region' global attribute matches the current region
				if 'region' in ds.attrs and ds.attrs['region'] == region:
					# If it matches, add the ncfile to the list for this region
					matching_ncfiles_dict[region].append(ncfile)

				# Close the dataset
				ds.close()

	return matching_ncfiles_dict, grdfingerprint_dict

'''
scale_ncfiles_by_fingerprint takes matching_ncfiles_dict and grdfingerprint_dict.
Loop over all the ice sources in grdfingerprint_dict. For each ice source,
loop over each region and evaluate whether there is a non-empty list of files
in matching_ncfiles_dict. If there is a matching file, load the fingerprint
for that region, and scale the matching netcdf files by the fingerprint. Add
the resulting scaled netcdf file to a dictionary that has keys of ice sources
and values that are summed arrays of scaled netcdf files. Return this dictionary.
'''

def scale_ncfiles_by_fingerprint(matching_ncfiles_dict, grdfingerprint_dict,
								 site_lats, site_lons, chunksize):
	
	nsites = len(site_lats)

	# Initialize the dictionary to store the scaled and summed netcdf files
	scaled_ncfiles_dict = {}

	# Loop over all the ice sources in grdfingerprint_dict
	for ice_source, regions in grdfingerprint_dict.items():
		# Initialize an empty array for this ice source
		scaled_ncfiles_dict[ice_source] = 0

		# Loop over each region
		for region, fingerprint in regions.items():
			# Check if there is a non-empty list of files in matching_ncfiles_dict
			if region in matching_ncfiles_dict and matching_ncfiles_dict[region]:
				# Load the fingerprint for that region
				regionfp = da.from_array(AssignFP(fingerprint, site_lats, site_lons), chunks=chunksize)

				# Loop over the matching netcdf files
				for ncfile in matching_ncfiles_dict[region]:

					# Load the netcdf file
					ds = xr.open_dataset(ncfile)

					# Initialize variable to hold the localized projections
					(nsamps, nregions, nyears) = ds['sea_level_change'].shape
					local_sl = da.zeros((nsamps, nyears, nsites), chunks=(-1,-1,chunksize))

					# Scale the netcdf file by the fingerprint
					local_sl += np.multiply.outer(ds['sea_level_change'], regionfp)

					# Add the resulting scaled netcdf file to the array for this ice source
					scaled_ncfiles_dict[ice_source] += local_sl
					targyears = ds['years']

					# Close the dataset
					ds.close()

	return scaled_ncfiles_dict, targyears


def write_localized_projections(scaled_ncfiles, site_ids, site_lats, site_lons, pipeline_id,
								scenario, baseyear, ice_source):

	# Define the missing value for the netCDF files
	nc_missing_value = np.nan #np.iinfo(np.int16).min

	# Create the xarray data structures for the localized projections
	ncvar_attributes = {"description": "Local SLR contributions according to emulandice2",
			"history": "Created " + time.ctime(time.time()),
			"source": "SLR Framework: emulandice2 workflow",
			"scenario": scenario,
			"baseyear": baseyear,
			"ice_source": ice_source}

	(nsamps, nregions, nyears) = scaled_ncfiles.shape
	rsl_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), scaled_ncfiles, 
											{"units":"mm", "missing_value":nc_missing_value}),
							"lat": (("locations"), site_lats),
							"lon": (("locations"), site_lons)},
		coords={"years": targyears, "locations": site_ids, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)

	rsl_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the post-processing stage for the emulandice glaciers SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
	parser.add_argument('--chunksize', help="Number of locations to process at a time [default=50]", type=int, default=50)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--grdfingerprintfile', help='YAML file that contains the fingerprints for each region')
	parser.add_argument('--scenario', help='Name of the scenario')
	parser.add_argument('--baseyear', help='Base year for the projections')
	parser.add_argument('--ncfiles', nargs='+', help='List of netCDF files')
	# Parse the arguments
	args = parser.parse_args()

	# Run the postprocessing for the parameters specified from the command line argument
	emulandice_postprocess(args.locationfile, args.chunksize, args.pipeline_id, args.ncfiles, args.grdfingerprintfile, args.scenario, args.baseyear)

	# Done
	exit()
