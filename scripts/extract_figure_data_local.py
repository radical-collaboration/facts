import numpy as np
import os
import sys
import re
import fnmatch
import argparse
from netCDF4 import Dataset

def extract_data(infilename, indir, outdir):

	# Generate this segments file name
	outfilename = re.sub(".nc", "_figuredata.nc", infilename)

	# Quantiles to extract from the data
	q = np.array([0.05, 0.17, 0.5, 0.83, 0.95])

	# Years of interest
	years = np.arange(2020,2301,10)

	# Open the netCDF file
	infile = os.path.join(indir, infilename)
	nc = Dataset(infile, 'r')

	# Get the quantiles
	qvar = nc.variables['quantiles'][:]
	ncyears = nc.variables['years'][:]

	# Find the indices that match up to the requested quantiles
	ncq_idx = np.flatnonzero(np.isin(qvar, q))
	ncq_idx = ncq_idx[np.arange(len(ncq_idx)-1)]  # Avoids spurious 95th percentile

	# Find the indices that match up to the requested years
	(_, years_idx, ncyears_idx) = np.intersect1d(years, ncyears, return_indices=True)

	# Get the data
	localsl_q = nc.variables['localSL_quantiles'][ncq_idx, :, ncyears_idx]
	ids = nc.variables['id'][:]
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
	site_dim = rootgrp.createDimension("nsites", localsl_q.shape[1])
	year_dim = rootgrp.createDimension("years", len(years_idx))
	q_dim = rootgrp.createDimension("quantiles", len(ncq_idx))

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
	rootgrp.description = nc_description
	rootgrp.history = nc_history
	rootgrp.source = nc_source
	lat_var.units = "Degrees North"
	lon_var.units = "Degrees East"
	localslq.units = "mm"

	# Put the data into the netcdf variables
	lat_var[:] = lats
	lon_var[:] = lons
	id_var[:] = ids
	year_var[:] = years[years_idx]
	q_var[:] = qvar[ncq_idx]
	localslq[:,:,:] = localsl_q

	# Close the netcdf
	rootgrp.close()


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Produce the figure data files from the FACTS localized projection files")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--indir', help="Input directory containing localized files (including total files)", required=True)
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
