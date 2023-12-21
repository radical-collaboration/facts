import sys
import argparse
import subprocess
import os

def emulandice_project(pipeline_id, ice_source, region, emu_name, emu_file, climate_data_file, scenario, nsamps, baseyear, 
					   seed, pyear_start, pyear_end, pyear_step):

	# Run the module using the FACTS forcing data
	arguments = [ice_source, region, emu_name, emu_file, climate_data_file, scenario, './', str(seed), pipeline_id]
	print(arguments)
	subprocess.run(["bash", "emulandice_steer.sh", *arguments])

	## MAYBE ABLE TO DO THIS WITHIN RPY2 WITH SOMETHING LIKE:

	# from rpy2.robjects import r
	# from rpy2.robjects.packages import importr

	# # Init the packrat environment
	# r.setwd('/path/to/packrat/environment') ## This should be the task directory
	# r.source('.Rprofile')

	# # use the packages it contains
	# importr('emulandice2')    
	# result = r.main()

	# BUT PERHAPS THEN NEED TO CALL emulandice_environment.sh before calling emulandice_project.py
	# THIS WOULD HAVE THE ADVANTAGE THAT r.main() could pass output as a return value rather
	# than having to write to a file and then read it back in if we want to do checks on it in this script

	## NEED TO MODIFY CHECKS SO THAT WORK WITH NETCDF OUTPUT

	# Get the output from the emulandice run
	# samples, targyears = ExtractProjections(emulandice_file)

	# Make sure we get the number of samples we expected
	#if nsamps != samples.shape[1]:
	#	raise Exception("Number of SLC projections does not match number of temperature trajectories: {} != {}".format(samples.shape[1], nsamps))


	return(None)

if __name__ == "__main__":

	# Initialize the argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the emulandice SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Add arguments for the resource and experiment configuration files
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--ice_source', help="Ice source: GIS, AIS or GLA", default='AIS', choices=['AIS','GIS','GLA'])
	parser.add_argument('--region', help="Ice source region: ALL for GIS/AIS and RGI01-RGI19 for GLA", default='ALL')
	parser.add_argument('--emu_name', help="Emulator name")
	parser.add_argument('--emu_file', help="Emulator file")
	parser.add_argument('--scenario', help="SSP Emissions scenario", default='ssp245')
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
	parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing
	emulandice_project(args.pipeline_id, args.ice_source, args.region, args.emu_name, args.emu_file, args.climate_data_file, 
					   args.scenario, args.nsamps, args.baseyear, 
					   args.seed, args.pyear_start, args.pyear_end, args.pyear_step)

	# Done
	sys.exit()
