import numpy as np
import os
import sys
import re
import argparse
import glob
#import time
from pathlib import Path
import xarray as xr
import dask.array as da


def FindFiles(indir):

	pattern = os.path.join(indir, "**", "*.nc")
	files = glob.iglob(pattern, recursive=True)

	return(files)


def GenerateRateFiles(indir, tempdir, chunksize, outdir=None):

	# Get a list of the files to convert
	infiles = FindFiles(indir)

	# Loop over the input files
	for this_file in infiles:

		# Generate the output file name
		if outdir is None:
			outdir = indir
		inpath = str(Path(this_file).parent)
		outpath = Path(re.sub(indir, outdir, inpath))
		Path.mkdir(outpath, parents=True, exist_ok=True)
		outfile = re.sub(indir, outdir, this_file)
		outfile = re.sub(r"\.nc", "_rates.nc", outfile)

		# Generate the rate file for this input file
		print(str(Path(this_file).name))
		#time.sleep(1)
		GenerateRateFile_dev(this_file, tempdir, chunksize, outfile)

	return(None)


'''
def GenerateRateFile(infile, chunksize, outfile):

	# Open the file
	#with xr.open_dataset(infile, chunks={'locations':chunksize}) as nc:
	with xr.open_dataset(infile) as nc:

		x_out = nc.differentiate("years")
		x_out = x_out.rename({"sea_level_change":"sea_level_change_rate"})
		x_out.sea_level_change_rate.attrs["units"] = "mm per year"

	# Missing value
	nc_missing_value = np.iinfo(np.int16).min

	# Write the rates to the output file
	x_out.to_netcdf(outfile, encoding={"sea_level_change_rate": {"scale_factor":0.1, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(None)
'''




# For merging multiple data sets with non-monotonically increasing/decreasing values along concat dimension
def open_mfdataset_merge_only(paths, **kwargs):
	paths = sorted(glob.iglob(paths))
	#return xr.merge([xr.open_dataset(path, **kwargs) for path in paths])
	return xr.concat([xr.open_dataset(path, **kwargs) for path in paths], dim="locations")


def GenerateRateFile_dev(infile, tempdir, chunksize, outfile):


	# Missing value
	nc_missing_value = np.iinfo(np.int16).min

	# Open the file
	with xr.open_dataset(infile) as nc:

		# Get the locations to loop over
		x_ids = nc['locations'].values
		nsites = len(x_ids)

		# Loop over the chunks of samples
		for i in np.arange(0,nsites,chunksize):

			# This chunk's temporary file name
			temp_filename = os.path.join(tempdir, "temp_rate_file_{0:05d}.nc".format(int(i/chunksize)))

			# If this temporary file exists, move on to the next chunk
			#if os.path.isfile(temp_filename):
			#	print("{} found, skipping to next chunk".format(temp_filename))
			#	continue

			# Get the sample indices
			top_idx = np.amin([i+chunksize, nsites])
			site_idx = np.arange(i,top_idx)

			# Subset the dataset for this chunk
			temp_out = nc.isel({'locations': site_idx})
			temp_out = temp_out.differentiate("years")
			temp_out = temp_out.rename({"sea_level_change":"sea_level_change_rate"})
			temp_out.sea_level_change_rate.attrs["units"] = "mm per year"

			# Write the rates to the output file
			temp_out.to_netcdf(temp_filename, encoding={"sea_level_change_rate": {"scale_factor":0.1, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Open the temporary data sets
	globstring = os.path.join(tempdir, "temp_rate_file_*.nc")
	combined = open_mfdataset_merge_only(globstring, chunks={"locations":chunksize})

	# Write the combined data out to the final netcdf file
	combined.to_netcdf(outfile, encoding={"sea_level_change_rate": {"scale_factor":0.1, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Remove the temporary files
	for i in np.arange(0,nsites,chunksize):
		os.remove(os.path.join(tempdir, "temp_rate_file_{0:05d}.nc".format(int(i/chunksize))))


	return(None)






if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Total up the contributors for a particular workflow.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--indir', help="Top level directory to search for input files", required=True)
	parser.add_argument('--outdir', help="Output directory for rate files [default=indir]", default=None)
	parser.add_argument('--tempdir', help="Temporary directory to hold intermediate rate files", required=True)
	parser.add_argument('--chunksize', help="Number of locations per chunk", default=50, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Total up the workflow in the provided directory
	GenerateRateFiles(args.indir, args.tempdir, args.chunksize, args.outdir)

	# Done
	sys.exit()
