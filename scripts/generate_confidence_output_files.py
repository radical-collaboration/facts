import numpy as np
import os
import sys
import re
import argparse
import fnmatch
import time
from pathlib import Path
import xarray as xr
import dask.array as da



def GetScenarios(dir):

	# Get the scenario names from the available scenario directories
	# Ignore any hidden directories (i.e. .DS_Store)
	pb1e_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_1e")) if not re.search(r"^\.", x)]
	pb1f_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_1f")) if not re.search(r"^\.", x)]
	pb2e_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_2e")) if not re.search(r"^\.", x)]
	pb2f_scenarios = [x for x in os.listdir(os.path.join(dir, "pb_2f")) if not re.search(r"^\.", x)]

	# Find the overlapping scenarios
	med_scenarios = list(set(pb1e_scenarios) & set(pb1f_scenarios))
	low_scenarios = list(set(pb2e_scenarios) & set(pb2f_scenarios))

	# Return the overlapping scenarios
	med_scenarios.sort()
	low_scenarios.sort()
	return(med_scenarios, low_scenarios)



def GetFiles(dir):

	# There should be a file for glaciers, landwaterstorage, oceandynamics,
	# AIS, GIS, and total...and for regional projections verticallandmotion.
	file_keys = ["glaciers", "landwaterstorage", "oceandynamics", "AIS", "GIS", "total", "verticallandmotion"]

	# Initialize list of matched file keys
	match_files = {}

	# Loop over the keys and find the associated files
	for this_key in file_keys:

		# Locate this file in the directory
		pattern = "*{}*.nc".format(this_key)
		this_file = fnmatch.filter(os.listdir(dir), pattern)

		# There should be only one match
		if len(this_file) == 1:

			# Append the match
			match_files[this_key] = os.path.join(dir, this_file[0])

		elif len(this_file) > 1:
			raise Exception("More than one file matched in {} for key {}".format(dir, this_key))

		else:

			match_files[this_key] = None

	# Return the dictionary of files
	return(match_files)



def MakeConfidenceFile(infile_e=None, infile_f=None, f_years=np.arange(2020,2301,10), outfile=None, is_rates=False, chunksize=50):

	# If both infile_e and infile_f are None, then there's no data for this component key.
	# Return and let the code move onto the next component key
	if infile_f is None and infile_e is None:
		return(1)

	# Variable names and attributes
	if is_rates:
		varname = "sea_level_change_rate"
		varscale = 0.1
	else:
		varname = "sea_level_change"
		varscale = 1.0

	# Open and subset the f file
	with xr.open_dataset(infile_f, chunks={"locations":chunksize}) as nc_f:
	#with xr.open_dataset(infile_f) as nc_f:
		nc_out = nc_f.sel(years=f_years)

	# Add the f file to the source list
	source_files = [infile_f]

	# If there's an e file, overlap it with the f file
	if infile_e is not None:
		with xr.open_dataset(infile_e, chunks={"locations":chunksize}) as nc_e:
		#with xr.open_dataset(infile_e) as nc_e:
			nc_out = nc_e.combine_first(nc_f.sel(years=f_years))

		# Append the e file to the source file list
		source_files.append(infile_e)

	# Define the missing value for the netCDF files
	nc_missing_value = np.iinfo(np.int16).min

	# Attributes for the output file
	nc_attrs = {"description": "Combined confidence output file for AR6 sea-level change projections",
			"history": "Created " + time.ctime(time.time()),
			"source": "Files Combined: {}".format(",".join(source_files))}

	# Put the attributes onto the output file
	nc_out.attrs = nc_attrs

	# Write the output file
	nc_out.to_netcdf(outfile, encoding={varname: {"scale_factor": varscale, "dtype": "i2", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	# Done
	return(None)


def GenerateConfidenceFiles(pboxdir, outdir, chunksize):

	# Are we working with values or rates?
	is_rates = True if re.search(r"rates", pboxdir) is not None else False

	# Get the overlapping scenarios for each confidence level
	med_scenarios, low_scenarios = GetScenarios(pboxdir)

	# If these are rate pboxes...
	if is_rates:

		# Medium-confidence: pb_1f through 2150
		# Low-confidence: pb_2f through 2300

		# Loop over the medium scenarios
		for this_scenario in med_scenarios:

			# Get the list of files for this scenario
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_1f", this_scenario))

			# Loop over the available components
			for this_key in pb1f_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "medium_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_medium_confidence_rates.nc".format(this_key, this_scenario))

				# Make the output file
				MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile, is_rates=is_rates, chunksize=chunksize)

		# Loop over the low scenarios
		for this_scenario in low_scenarios:

			# Get the list of files for this scenario
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_2f", this_scenario))

			# Loop over the available components
			for this_key in pb1f_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "low_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_low_confidence_rates.nc".format(this_key, this_scenario))

				# Make the output file
				MakeConfidenceFile(infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2301,10), outfile=outfile, is_rates=is_rates, chunksize=chunksize)


	# These are value files...
	else:

		# Medium-confidence: pb_1e through 2100, pb_1f through 2150
		# Low-confidence: pb_2e through 2100, pb_2f through 2300

		# Loop over the medium scenarios
		for this_scenario in med_scenarios:

			# Get the list of files for this scenario
			pb1e_infiles = GetFiles(os.path.join(pboxdir, "pb_1e", this_scenario))
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_1f", this_scenario))

			# Loop over the available components
			for this_key in pb1e_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "medium_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_medium_confidence_values.nc".format(this_key, this_scenario))

				# Make the output file
				MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2151,10), outfile=outfile, chunksize=chunksize)

		# Loop over the low scenarios
		for this_scenario in low_scenarios:

			# Get the list of files for this scenario
			pb1e_infiles = GetFiles(os.path.join(pboxdir, "pb_2e", this_scenario))
			pb1f_infiles = GetFiles(os.path.join(pboxdir, "pb_2f", this_scenario))

			# Loop over the available components
			for this_key in pb1e_infiles.keys():

				# Define the output file name
				outpath = Path(os.path.join(outdir, "low_confidence", this_scenario))
				Path.mkdir(outpath, parents=True, exist_ok=True)
				outfile = os.path.join(outpath, "{}_{}_low_confidence_values.nc".format(this_key, this_scenario))

				# Make the output file
				MakeConfidenceFile(infile_e=pb1e_infiles[this_key], infile_f=pb1f_infiles[this_key], f_years=np.arange(2020,2301,10), outfile=outfile, chunksize=chunksize)

	# Done
	return(None)




if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Generate the medium- and low-confidence output files for AR6 projections.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pboxdir', help="Directory that contains the pbox workflows (pb_1e, pb_1f, pb_2e, and pb_2f)", required=True)
	parser.add_argument('--outdir', help="Top-level output directory for output files", required=True)
	parser.add_argument('--chunksize', help="Number of locations per chunk", default=50, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Total up the workflow in the provided directory
	GenerateConfidenceFiles(args.pboxdir, args.outdir, args.chunksize)

	# Done
	sys.exit()
