import numpy as np
import os
import sys
import re
import fnmatch
import argparse
from netCDF4 import Dataset

def extract_data(infilename, indir, outdir):

	# Initialize valid variable names and units
	varname = None
	varunit = None
	valid_varnames = ["sea_level_change", "sea_level_change_rate"]
	valid_varunits = {"sea_level_change": "mm", "sea_level_change_rate": "mm per year"}
	valid_varscale = {"sea_level_change": 1.0, "sea_level_change_rate": 0.1}

	# Generate this segments file name
	outfilename = re.sub(".nc", "_figuredata.nc", infilename)

	# Quantiles to extract from the data
	q = np.array([0.05, 0.17, 0.5, 0.83, 0.95])

	# Years of interest
	years = np.arange(2020,2301,10)

	# Open the netCDF file
	infile = os.path.join(indir, infilename)
	nc = Dataset(infile, 'r')

	# Check for a valid variable. If none exists, print message to STDOUT and
	# return to main
	for this_varname in valid_varnames:
		if this_varname in nc.variables.keys():
			varname = this_varname
			varunit = valid_varunits[this_varname]
			varscale = valid_varscale[this_varname]
			break
	if varname is None:
		nc.close()
		print("No valid variable name exists in {}".format(infile))
		return(1)

	# Get the quantiles
	qvar = nc.variables['quantiles'][:]
	ncyears = nc.variables['years'][:]

	# Find the indices that match up to the requested quantiles
	ncq_idx = np.flatnonzero(np.isin(qvar, q))

	# Find the indices that match up to the requested years
	(_, years_idx, ncyears_idx) = np.intersect1d(years, ncyears, return_indices=True)

	# Get the data
	localsl_q = nc.variables[varname][ncq_idx, ncyears_idx, :]
	ids = nc.variables['locations'][:]
	lats = nc.variables['lat'][:]
	lons = nc.variables['lon'][:]

	nc_description = nc.__dict__['description']
	nc_history = nc.__dict__['history']
	nc_source = nc.__dict__['source']

	# Close the input netcdf file
	nc.close()

	# Define the output file
	outfile = os.path.join(outdir, outfilename)

	# Write the combined projections to a netcdf file
	rootgrp = Dataset(outfile, "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("locations", localsl_q.shape[2])
	year_dim = rootgrp.createDimension("years", len(years_idx))
	q_dim = rootgrp.createDimension("quantiles", len(ncq_idx))

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))
	id_var = rootgrp.createVariable("locations", "i4", ("locations",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable(varname, "i2", ("quantiles", "years", "locations"), zlib=True, complevel=4)
	localslq.missing_value = np.iinfo(np.int16).min
	localslq.scale_factor = varscale

	# Assign attributes
	rootgrp.description = nc_description
	rootgrp.history = nc_history
	rootgrp.source = nc_source
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	localslq.units = varunit

	# Put the data into the netcdf variables
	lat_var[:] = lats
	lon_var[:] = lons
	id_var[:] = ids
	year_var[:] = years[years_idx]
	q_var[:] = qvar[ncq_idx]
	localslq[:,:,:] = localsl_q

	# Close the netcdf
	rootgrp.close()

	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Produce the figure data files from the FACTS projection files")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--indir', help="Input directory containing projection files (including total files)", required=True)
	parser.add_argument('--outdir', help="Output directory to contain the figure data", required=True)

	# Parse the arguments
	args = parser.parse_args()

	# Loop through the files
	for this_file in os.listdir(args.indir):

		if not re.search(".nc", this_file):
			continue

		# Extract the data from the file
		extract_data(this_file, args.indir, args.outdir)

	# Done
	sys.exit()
