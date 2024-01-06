import sys
import argparse
import subprocess
import os
import time
import xarray as xr
import numpy as np
import dask
import scipy

def emulandice_project(pipeline_id, ice_source, regions, emu_file, climate_data_file, scenario, baseyear, 
					   seed, pyear_start, pyear_end, pyear_step):

	# Run the module using the FACTS forcing data
	if len(regions) != len(emu_file):
		raise Exception("Number of regions and emulator files must be the same")

	targyears = np.arange(pyear_start,pyear_end+.5,pyear_step)
	lazy_regions = [dask.delayed(RebaseSamples(emulandice_steer(ice_source, regions[ii], emu_file[ii],
															  climate_data_file, scenario, './', seed, pipeline_id),
															  targyears, baseyear)) for ii in np.arange(len(regions))] 
	run_regions = dask.compute(*lazy_regions)

	if len(regions) > 1:
		outfile = pipeline_id + "_ALL_globalsl.nc"
		# does all-region outfile already exist? if not, produce it
		if not os.path.exists(outfile):
			TotalSamples(run_regions, outfile, 50, ice_source)

	return(run_regions)

def RebaseSamples(ncfile,targyears,baseyear):
	print('Rebasing ' + ncfile + '...')
	ds = xr.open_dataset(ncfile)
	attrs=ds.attrs
	ds = ds.interp(years=targyears)-ds.interp(years=baseyear)
	attrs['baseyear'] = baseyear
	ds.attrs = attrs
	print(ds.attrs)
	ds.to_netcdf(ncfile,encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4}})
	return(ncfile)

def emulandice_steer(ice_source, region, emu_file, climate_data_file, scenario, outdir, seed, pipeline_id ):
	arguments = [ice_source, region, emu_file, climate_data_file, scenario, outdir, str(seed), pipeline_id]
	subprocess.run(["bash", "emulandice_steer.sh", *arguments])
	outfile = pipeline_id + "_" + region + "_globalsl.nc"
	return(outfile)

def TotalSamples(infiles, outfile, chunksize, ice_source):
    # Reads in multiple files (delayed) and tries combine along
    # common dimensions and a new "file" dimension.
	ds = xr.open_mfdataset(
		infiles, 
	    combine="nested", 
	    concat_dim="file", 
	    chunks={"locations":chunksize},
	)

	# Sums everything across the new "file" dimension.
	total_out = ds[["sea_level_change"]].sum(dim="file")
	# Add "lat" and "lon" as data variable in output, pulling values from the first file.
	total_out["lat"] = ds["lat"].isel(file=0)
	total_out["lon"] = ds["lon"].isel(file=0)

	# Attributes for the total file
	total_out.attrs = {
		"description": "Total " + ice_source + " sea-level change",
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
	
	total_out.to_netcdf(outfile, encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})

	return(outfile)


if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the emulandice SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--ice_source', help="Ice source: GIS, AIS or GLA", default='AIS', choices=['AIS','GIS','GLA'])
	parser.add_argument('--region', nargs='+', help="Ice source region: ALL for GIS/AIS and RGI01-RGI19 for GLA", default='ALL')
	parser.add_argument('--emu_file', nargs='+', help="Emulator file")
	parser.add_argument('--scenario', help="SSP Emissions scenario", default='ssp245')
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_project(args.pipeline_id, args.ice_source, args.region, args.emu_file, args.climate_data_file, 
					   args.scenario, args.baseyear, 
					   args.seed, args.pyear_start, args.pyear_end, args.pyear_step)

	# Done
	sys.exit()
