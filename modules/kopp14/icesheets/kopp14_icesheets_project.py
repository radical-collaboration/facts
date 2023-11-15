import pickle
import sys
import os
import argparse
import numpy as np
from SampleISDists import SampleISDists
from ProjectGSL import ProjectGSL
import cholcov
import time
from netCDF4 import Dataset

''' kopp14_project_icesheets.py

This script runs the ice sheet projection task for the Kopp 2014 workflow. This task
requires data produced from the pre-processing and fit/calibration stages. This task
also requires at the time of call the number of samples to produce and the names of the
input files generated from the previous stages.

Parameters:
nsamps = Number of samples to produce
seed = The seed to use for the random number gerneration
pipeline_id = Unique identifier for the pipeline running this code

Output:
- Pickled file containing "nsamps" of global sea-level rise due to ice sheets (internal)
- NetCDF files for GIS, EAIS, and WAIS global contributions (external)

'''

def kopp14_project_icesheets(nsamps, seed, baseyear, pyear_start, pyear_end, pyear_step, pipeline_id):

	## Read in the fitted parameters from parfile
	# Open the file
	parfile = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(parfile, 'rb')
	except:
		print("Cannot open parfile\n")

	# Extract the data from the file
	my_data = pickle.load(f)
	batheteais = my_data['batheteais']
	bathetwais = my_data['bathetwais']
	bathetgis = my_data['bathetgis']
	arthetais = my_data['arthetais']
	arthetgis = my_data['arthetgis']
	islastdecade = my_data['islastdecade']
	f.close()

	## Read in the correlation information
	corfile = "{}_corr.pkl".format(pipeline_id)
	try:
		f = open(corfile, 'rb')
	except:
		print("Cannot open corfile\n")

	# Extract the data
	my_data = pickle.load(f)
	bacorris = my_data['bacorris']
	arcorris = my_data['arcorris']
	f.close()

	# Generate samples of ice sheet accelerations
	sigmas = np.array([bathetgis[2], bathetwais[2], batheteais[2]])
	mus = np.array([bathetgis[1], bathetwais[1], batheteais[1]])
	offsets = np.array([bathetgis[0], bathetwais[0], batheteais[0]])
	baissamps = SampleISDists(nsamps, sigmas, mus, offsets, islastdecade, bacorris, seed)

	sigmas = np.array([arthetgis[2], arthetais[2]])
	mus = np.array([arthetgis[1], arthetais[1]])
	offsets = np.array([arthetgis[0], arthetais[0]])
	arislastdecade = np.array([islastdecade[0], islastdecade[1]+islastdecade[2]])
	arissamps = SampleISDists(nsamps, sigmas, mus, offsets, arislastdecade, arcorris, seed+1234)

	# Project global sea-level rise over time
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	targyears = np.union1d(targyears, baseyear)
	(arsamps, basamps, hysamps) = ProjectGSL(baissamps, arissamps, islastdecade, targyears)

	# Reference these projections to the base year
	baseyear_idx = np.flatnonzero(targyears == baseyear)
	arsamps = arsamps - arsamps[:,baseyear_idx,:]
	basamps = basamps - basamps[:,baseyear_idx,:]
	hysamps = hysamps - hysamps[:,baseyear_idx,:]

	# Put the results into a dictionary
	output = {'arsamps': arsamps, 'basamps': basamps, 'hysamps': hysamps, 'targyears': targyears, 'baseyear': baseyear}

	# Write the results to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the global projection output files
	writeNetCDF(hysamps[:,:,0], pipeline_id, "GIS", targyears, baseyear, nsamps)
	writeNetCDF(hysamps[:,:,1], pipeline_id, "WAIS", targyears, baseyear, nsamps)
	writeNetCDF(hysamps[:,:,2], pipeline_id, "EAIS", targyears, baseyear, nsamps)
	writeNetCDF(hysamps[:,:,1]+hysamps[:,:,2], pipeline_id, "AIS", targyears, baseyear, nsamps)


	# Done
	return(0)



def writeNetCDF(data, pipeline_id, icesheet_name, targyears, baseyear, nsamps):

	# Write the localized projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_{}_globalsl.nc".format(pipeline_id, icesheet_name))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = "Global SLR contribution from {} according to Kopp 2014 workflow".format(icesheet_name)
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}, Baseyear: {1}".format(pipeline_id, baseyear)
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nsamps)
	samps[:,:,:] = data[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()

	# Done
	return(0)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--baseyear', help="Base or reference year for projetions [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Make sure the base year and target years are within data limits for this module
	if(args.baseyear < 2000):
		raise Exception("Base year cannot be less than year 2000: baseyear = {}".format(args.baseyear))
	if(args.baseyear > 2300):
		raise Exception("Base year cannot be greater than year 2300: baseyear = {}".format(args.baseyear))
	if(args.pyear_start < 2000):
		raise Exception("Projection year cannot be less than year 2000: pyear_start = {}".format(args.pyear_start))
	if(args.pyear_end > 2300):
		raise Exception("Projection year cannot be greater than year 2300: pyear_end = {}".format(args.pyear_end))

	# Make sure the target year stepping is positive
	if(args.pyear_step < 1):
		raise Exception("Projection year step must be greater than 0: pyear_step = {}".format(args.pyear_step))

	# Run the projection process for the parameters specified from the command line argument
	kopp14_project_icesheets(args.nsamps, args.seed, args.baseyear, args.pyear_start, args.pyear_end, args.pyear_step, args.pipeline_id)

	# Done
	exit()
