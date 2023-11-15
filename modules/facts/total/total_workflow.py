import numpy as np
import os
import re
import time
import argparse
import shutil
import yaml
import xarray as xr
import dask.array as da


def TotalSamplesInDirectory(directory, pyear_start, pyear_end, pyear_step, chunksize):

	# Output directory
	outdir = os.path.dirname(__file__)

	# Define the output file
	outfilename = "total.workflow.nc"
	outfile = os.path.join(outdir, outfilename)

	# Get the list of input files
	infiles0 = [file for file in os.listdir(directory) if file.endswith(".nc")]
	
	# Append the input directory to the input file name
	infiles = []
	for infile in infiles0:
		infiles.append(os.path.join(directory, infile))

	# Define the years of interest
	targyears = xr.DataArray(np.arange(pyear_start, pyear_end+1, pyear_step), dims="years")

	r=TotalSamples(infiles, outfile, targyears, chunksize)

	# Put a copy of the total file back into the shared directory
	shutil.copy2(outfile, directory)

	return(r)


def TotalSampleInWorkflow(wfcfg, directory, targyears, workflow, scale, chunksize=50, experiment_name=None):
		# Define the output file
		outdir = os.path.dirname(__file__)
		outfilename = "total.workflow." + workflow + "." + scale + ".nc"
		if len(experiment_name)>0:
			outfilename = experiment_name + "." + outfilename
		outfile = os.path.join(outdir, outfilename)

		# Get the list of input files
		infiles = []
		r=[]
		for this_file in wfcfg[workflow][scale]:
			infiles.append(os.path.join(directory,this_file))

		if len(infiles)>0:
			rout=TotalSamples(infiles, outfile, targyears, chunksize)
			r.append(rout)

			# Put a copy of the total file back into the shared directory
			shutil.copy2(outfile, directory)
		return(r)


def TotalSamplesInWorkflows(directory, pyear_start, pyear_end, pyear_step, chunksize, wfcfg_file="workflows.yml", workflow=[""], scale=[""], experiment_name=None):

	# Define the years of interest
	targyears0 = xr.DataArray(np.arange(pyear_start, pyear_end+1, pyear_step), dims="years")

	# read yaml file
	with open(wfcfg_file, 'r') as fp:
		wfcfg = yaml.safe_load(fp)

	# loop through workflows and scopes
	r=[]
	for this_workflow in wfcfg:

		targyears = targyears0

		if 'options' in wfcfg[this_workflow].keys():
			if 'pyear_end' in wfcfg[this_workflow]['options'].keys():
				targyears = xr.DataArray(np.arange(pyear_start, min([pyear_end,wfcfg[this_workflow]['options']['pyear_end']])+1, pyear_step), dims="years")

		for this_scale in wfcfg[this_workflow]:
			if this_scale in {'options'}:
				continue
			
			if len(workflow[0])>0:
				if this_workflow in workflow:
					if len(scale[0])>0:
						if this_scale in scale:
							r.append(TotalSampleInWorkflow(wfcfg, directory, targyears, this_workflow, this_scale, chunksize, experiment_name=experiment_name))
						else:
							r.append(TotalSampleInWorkflow(wfcfg, directory, targyears, this_workflow, this_scale, chunksize, experiment_name=experiment_name))
					else:
						r.append(TotalSampleInWorkflow(wfcfg, directory, targyears, this_workflow, this_scale, chunksize, experiment_name=experiment_name))
			else:
				r.append(TotalSampleInWorkflow(wfcfg, directory, targyears, this_workflow, this_scale, chunksize, experiment_name=experiment_name))

	return(r)

def TotalSamples(infiles, outfile, targyears, chunksize):
	# Skip this file if it appears to be a total file already
	# SBM: FWIW, this might not work as intended because path is appended to file names in funcs above.
	target_infiles = [fl for fl in infiles if not re.search(r"^total", fl)]

    # Reads in multiple files (delayed) and tries combine along
    # common dimensions and a new "file" dimension.
	ds = xr.open_mfdataset(
		target_infiles, 
	    combine="nested", 
	    concat_dim="file", 
	    chunks={"locations":chunksize},
	)
	ds = ds.sel(years=targyears)
	# Sums everything across the new "file" dimension.
	total_out = ds[["sea_level_change"]].sum(dim="file")
	# Add "lat" and "lon" as data variable in output, pulling values from the first file.
	total_out["lat"] = ds["lat"].isel(file=0)
	total_out["lon"] = ds["lon"].isel(file=0)

	# Attributes for the total file
	total_out.attrs = {
		"description": "Total sea-level change for workflow",
		"history": "Created " + time.ctime(time.time()),
		"source": "FACTS: Post-processed total among available contributors: {}".format(",".join(infiles)),
	}

	# Define the missing value for the netCDF files
	nc_missing_value = np.nan #np.iinfo(np.int16).min
	total_out["sea_level_change"].attrs = {
		"units": "mm", 
		"missing_value": nc_missing_value
	}

	# Write the total to an output file.
    # This actually carries out the delayed calculations and operations.
    # SBM: FYI Double check the numbers to ensure everything is summing across dims correctly.
    # SBM: FYI Also, check to see if output as something huge like float64.
	total_out.to_netcdf(outfile, encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(outfile)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Total up the contributors for a particular workflow.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--directory', help="Directory containing files to aggregate", required=True)
	parser.add_argument('--workflows', help="Use a workflows.yml file listing files to aggregate", type=str, default=None)
	parser.add_argument('--workflow', help="Workflow to run", type=str, default="")
	parser.add_argument('--scale', help="Scale to run", type=str, default="")
	parser.add_argument('--experiment_name', help="Experiment Name", type=str, default="")
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--chunksize', help="Number of locations per chunk", default=50, type=int)

	# Parse the arguments
	args = parser.parse_args()

	if args.workflows or (len(args.workflow)>0) or (len(args.scale)>0):
		# Total up by workflow
		TotalSamplesInWorkflows(args.directory,  args.pyear_start, args.pyear_end, args.pyear_step, args.chunksize, wfcfg_file= args.workflows, workflow=[args.workflow], scale=[args.scale], experiment_name = args.experiment_name)
	else:
		# Total up the workflow in the provided directory
		TotalSamplesInDirectory(args.directory, args.pyear_start, args.pyear_end, args.pyear_step, args.chunksize)

	exit()
