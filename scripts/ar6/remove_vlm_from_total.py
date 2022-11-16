import os
import sys
import argparse
import re
import fnmatch
import numpy as np
import xarray as xr
import pickle
import shutil


def RemoveVLM(vlm_file, total_file, outfile, chunksize):

	# Initialize the variable name to extract
	varname = "sea_level_change"
	varscale = 1.0
	varunits = "mm"
	vardtype = np.int16
	nc_missing_value = np.iinfo(np.int16).min

	# Extract the VLM data
	with xr.open_dataset(vlm_file) as vlmnc:
		vlm_years = vlmnc.variables['years']

	# Extract the total data
	with xr.open_dataset(total_file) as totalnc:
		total_years = totalnc.variables['years']
		lats = totalnc['lat']
		lons = totalnc['lon']
		total_attrs = totalnc.attrs
		total_attrs['source'] = "FACTS - The Framework for Assessing Changes To Sea-level"

	# Check the overlap in years
	overlap_years = np.intersect1d(vlm_years, total_years)

	# Extract the VLM data
	with xr.open_dataset(vlm_file, chunks={'locations':chunksize}) as vlmnc:
		if not varname in vlmnc.keys():
			varname = "sea_level_change_rate"
			varscale = 0.1
			varunits = "mm per year"
			vardtype = np.float16
		vlm = vlmnc[varname].sel({'years': overlap_years}).astype(vardtype)

	# Extract the total data
	with xr.open_dataset(total_file, chunks={'locations':chunksize}) as totalnc:
		total = totalnc[varname].sel({'years': overlap_years}).astype(vardtype)

	# Subtract the two
	total_new = total - vlm

	# Append the units to the sea-level change variable
	total_new.attrs['units'] = varunits

	# Merge the coordinates with the data
	merged_data = xr.merge([total_new, lats, lons])
	merged_data.attrs = total_attrs

	# Write the output file
	merged_data.to_netcdf(outfile, encoding={varname: {"scale_factor": varscale, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)

def main(indir, outdir, chunksize):

	# Get the files from the input directory
	infiles = os.listdir(indir)
	infiles = fnmatch.filter(infiles, "*.nc")

	# Initialize
	vlm_file = None
	total_file = None

	# Loop over the input files
	for infile in infiles:

		# Note the vlm and total files
		if re.search("verticallandmotion", infile):
			vlm_file = infile
			continue
		if re.search("total", infile):
			total_file = infile
			continue

		# This input file is neither vlm or a total, so copy it over
		shutil.copy2(os.path.join(indir, infile), os.path.join(outdir, infile))

	# Remove the vlm from total if possible
	if vlm_file is not None and total_file is not None:
		RemoveVLM(os.path.join(indir, vlm_file), os.path.join(indir, total_file), os.path.join(outdir, total_file), chunksize)
	else:
		raise Exception("No VLM or Total file found")

	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Produces projections for a non-vlm workflow based on a workflow with vlm included.")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--indir', help="Directory that contains output from a FACTS projection workflow (with VLM)", required=True)
	parser.add_argument('--outdir', help="Output directory", required=True)
	parser.add_argument('--chunksize', help="Chunk size to use when reading in data", type=int, default=-1)

	# Parse the arguments
	args = parser.parse_args()


	main(args.indir, args.outdir, args.chunksize)

	sys.exit()
