import os
import sys
import argparse
import re
import time
import numpy as np
from netCDF4 import Dataset


def LoadInfiles(infiles, years):

	# Initialize the return variables
	localsl_q = []
	file_list = []
	ids = None
	lats = None
	lons = None
	qvar = None

	# Initialize the first file flag
	first_file = True

	# Loop over the input files
	for infile in infiles:

		# Open the input file for reading
		nc = Dataset(infile, "r")

		# Find which years in the projections match the requested years
		ncyears = nc.variables['years'][:]
		_, _, years_idx = np.intersect1d(years, ncyears, return_indices=True)

		# Extract the projection quantiles
		localsl_q.append(nc.variables['localSL_quantiles'][:,:,years_idx])

		# If this is the first file, extract the non-pbox information that gets passed
		# to the final pbox component file
		if first_file:
			ids = nc.variables['id'][:]
			lats = nc.variables['lat'][:]
			lons = nc.variables['lon'][:]
			qvar = nc.variables['quantiles'][:]

			# Set the first file flag to false
			first_file = False

		# Close the netcdf file
		nc.close()

	# Convert everything into numpy arrays
	localsl_q = np.array(localsl_q)
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	qvar = np.array(qvar)

	# Return everything
	return(localsl_q, ids, lats, lons, qvar)


def main(infiles, outfile, pyear_start=2020, pyear_end=2100, pyear_step=10):

	# Read in all the infiles and generate a large data array
	years = np.arange(pyear_start, pyear_end+1, pyear_step)
	component_data, ids, lats, lons, qvar = LoadInfiles(infiles, years)

	# Make the infile string for the output netcdf file
	infile_string = ", ".join([os.path.basename(x) for x in infiles])

	# Which indices are the median, above the median, and below the median
	median_idx = np.flatnonzero(qvar == 0.5)
	above_idx = np.arange(median_idx + 1, len(qvar))
	below_idx = np.arange(median_idx)

	# Initialize the output pbox array
	pbox = np.full(component_data.shape[1:], np.nan)

	# Use the mean of the medians of the components as the median for the pbox
	pbox[median_idx,:,:] = np.mean(component_data[:,median_idx,:,:], axis=0)

	# Use the minimum for below the median and maximum for above
	pbox[below_idx,:,:] = np.amin(component_data[:,below_idx,:,:], axis=0)
	pbox[above_idx,:,:] = np.amax(component_data[:,above_idx,:,:], axis=0)

	# Write the output netcdf file
	rootgrp = Dataset(outfile, "w", format="NETCDF4")

	# Define Dimensions
	site_dim = rootgrp.createDimension("nsites", pbox.shape[1])
	year_dim = rootgrp.createDimension("years", pbox.shape[2])
	q_dim = rootgrp.createDimension("quantiles", pbox.shape[0])

	# Populate dimension variables
	lat_var = rootgrp.createVariable("lat", "f4", ("nsites",))
	lon_var = rootgrp.createVariable("lon", "f4", ("nsites",))
	id_var = rootgrp.createVariable("id", "i4", ("nsites",))
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	q_var = rootgrp.createVariable("quantiles", "f4", ("quantiles",))

	# Create a data variable
	localslq = rootgrp.createVariable("localSL_quantiles", "i2", ("quantiles", "nsites", "years"), zlib=True, complevel=4)
	localslq.missing_value = np.iinfo(np.int16).min

	# Assign attributes
	rootgrp.description = "Pbox generated from a localized sea-level change projection workflow from FACTS"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "Generated with files: {}".format(infile_string)
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	localslq.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = lats
	lon_var[:] = lons
	id_var[:] = ids
	year_var[:] = years
	q_var[:] = qvar
	localslq[:,:,:] = pbox

	# Close the netcdf
	rootgrp.close()




if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Produce the pbox for the localized projection component files provided")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--infiles', help="Input files over which pbox is calculated", nargs="+", required=True)
	parser.add_argument('--outfile', help="Output file for the pbox", required=True)
	parser.add_argument('--pyear_start', help="Projection year start [default = 2020]", type=int, default=2020)
	parser.add_argument('--pyear_end', help="Projection year end [default = 2100]", type=int, default=2100)
	parser.add_argument('--pyear_step', help="Projection year step [default = 10]", type=int, default=10)

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	main(args.infiles, args.outfile, args.pyear_start, args.pyear_end, args.pyear_step)

	sys.exit()
