# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19

import os
import numpy as np
import pickle
import argparse
import time
import re
from netCDF4 import Dataset


class ProjectionError(Exception):
	pass


def project_glacier1(it,factor,exponent):
	
	# Return projection of glacier contribution by one glacier method
	scale=1e-3 # mm to m
	it[it < 0] = 0.0
	return scale*factor*(it**exponent)
	
def ar5_project_glaciersfair(rng_seed, nmsamps, ntsamps, pipeline_id):
	
	# Define the target years
	targyears = np.arange(2020,2101,10)
	
	# Load the preprocessed data
	data_file = "{}_data.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file\n")
	
	# Extract the data variables
	my_data = pickle.load(f)
	f.close()
	
	inttemp_mean = my_data['inttemp_mean']
	inttemp_sd = my_data['inttemp_sd']
	data_years = my_data['data_years']
	startyr = my_data["startyr"]
	scenario = my_data['scenario']
	#include_models = my_data['include_models']
	#include_scenarios = my_data['include_scenarios']
	#nmodels = len(include_models)
	
	# Load the fit data
	data_file = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open fit file\n")
	
	# Extract the fit variables
	my_fit = pickle.load(f)
	f.close()
	
	dmz = my_fit['dmz']
	#cvgl = my_fit['cvgl']
	glmass = my_fit['glmass']
	glparm = my_fit['glparm']
	
	# Produce the included model string
	#model_string_pieces = ["{0}-{1}".format(include_models[x], include_scenarios[x]) for x in np.arange(nmodels)]
	#model_string = "Models and scenarios included: " + ", ".join(model_string_pieces)
	model_string = "Temperature input from FAIR model ensemble"
	
	# Subset the data to the target years
	year_idx = np.isin(data_years, targyears)
	inttemp_mean = inttemp_mean[year_idx]
	inttemp_sd = inttemp_sd[year_idx]
	data_years = data_years[year_idx]
	
	# Set the seed for the random number generator
	np.random.seed(rng_seed)

	# Generate perfectly correlated samples
	z=np.random.standard_normal(ntsamps)[:,np.newaxis]
	
	# For each quantity, mean + standard deviation * normal random number
	zit=inttemp_mean + (inttemp_sd * z)

	# Number of realizations
	nr = nmsamps
	
	# Number of years in the data record
	nyr = len(data_years)
	
	# number of glacier methods
	ngl=len(glparm)
	
	# Stop if the number of realizations requested cannot evenly distributed over methods
	if nr%ngl:
		raise ProjectionError('number of realisations '+\
			'must be a multiple of number of glacier methods')
	
	# number of realisations per glacier method
	nrpergl=nr/ngl 
	
	# Generate samples for methodologies
	r=np.random.standard_normal(nr)
	
	# Initialize the data structure to hold the glacier samples
	glacier = np.full((nmsamps, ntsamps, nyr), np.nan)

	# Make an ensemble of projections for each method
	for igl in np.arange(ngl):
		
		# Load the methodological error for this method
		cvgl = glparm[igl]['cvgl']
		
		# glacier projection for this method using the mean temperature timeseries
		mgl=project_glacier1(inttemp_mean,glparm[igl]['factor'],glparm[igl]['exponent'])
		
		# glacier projections for this method with the ensemble of timeseries
		zgl=project_glacier1(zit,glparm[igl]['factor'],glparm[igl]['exponent'])
		
		# Store these samples
		ifirst=int(igl*nrpergl)
		ilast=int(ifirst+nrpergl)
		glacier[ifirst:ilast,...] = zgl[np.newaxis,:,:] + (mgl * r[ifirst:ilast,np.newaxis] * cvgl)[:,np.newaxis,:]
		
	
	# Add the contribution from the glaciers between the start year and
	# the AR5 reference period
	glacier += dmz
	glacier[glacier > glmass] = glmass
	
	# Flatten the sample data structure
	glacier = glacier.reshape(-1, glacier.shape[-1]) * 1000  # Convert to mm
	total_glac_samps = glacier.T

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(data_years))
	samp_dim = rootgrp.createDimension("samples", ntsamps*nmsamps)

	# Populate dimension variables
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	samps = rootgrp.createVariable("samps", "f4", ("years", "samples"), zlib=True, least_significant_digit=2)	
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from glaciers and ice caps according to AR5 method (calibrated to GMIP)"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}. ".format(pipeline_id, scenario) + model_string
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(0,nmsamps*ntsamps)
	samps[:,:] = total_glac_samps

	# Close the netcdf
	rootgrp.close()	
	
	# Load in the glacier fraction data-------------------------------------------
	# Note: There's were derived from the Kopp14 workflow.  Apparently, glacier region 4
	#		is not represented and region 7 is represented twice.  I'm not sure why, but
	#		it's consistent with the K14 workflow, so it's been carried over to this.
	
	# Initialize the data structures
	glac_frac = []
	glac_region_names = []
	
	# Open the glacier fraction file
	glac_frac_file = os.path.join(os.path.dirname(__file__), "glacier_fraction.txt")
	with open(glac_frac_file, 'r') as fp:
		
		# Get the fraction years from the header line
		header_items = re.split(",\s*", fp.readline())
		glac_frac_years = np.array([int(x) for x in header_items[1:]])
		
		# Read in the rest of the files
		for line in fp:
			line = line.rstrip()
			
			# Split the line into the region name and the fractions then append to data structures
			line_parts = re.split(",\s*", line)
			glac_region_names.append(line_parts[0])
			glac_frac.append([float(x) for x in line_parts[1:]])
	
	# Convert the fraction data structure into a numpy array
	glac_frac = np.array(glac_frac)
	
	# Subset the fraction data to the years of interest
	year_idx = np.isin(glac_frac_years, data_years)
	glac_frac = glac_frac[:,year_idx]
	
	# Reshape the samples and fraction data structures for broadcasting
	glacier = glacier[:,np.newaxis,:]
	glac_frac = glac_frac[np.newaxis,:,:]
	
	# Apply the regional fractions to the global projections
	gicsamps = glacier * glac_frac
	
	# Save the global glacier and ice caps projections to a pickle
	output = {"gicsamps": gicsamps, "glac_region_names": glac_region_names, "data_years": data_years}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(0)

	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glacier projection stage for the AR5 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nmsamps', help="Number of method samples to generate [default=200]", default=200, type=int)
	parser.add_argument('--ntsamps', help="Number of climate samples to generate [default=100]", default=100, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	ar5_project_glaciersfair(args.seed, args.nmsamps, args.ntsamps, args.pipeline_id)
	
	exit()
