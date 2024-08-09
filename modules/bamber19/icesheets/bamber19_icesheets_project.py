import numpy as np
import argparse
import pickle
import os
import re
import time
import h5py
from netCDF4 import Dataset

''' bamber19_project_icesheets.py

This is the projection stage for the Bamber et al. 2019 ice sheet component of the IPCC AR6 module set.

Parameters:
nsamps              Number of samples to produce
replace             Allow sampling with replacement
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline
climate_data_file   FAIR-outputted climate data file


Output:
"{pipeline_id}_projections.pkl" = User requested samples
"{pipeline_id}_(AIS|EAIS|WAIS|GIS)_globalsl.nc" = Sampled global projections in netCDF file

'''

def bamber19_project_icesheets(nsamps, pipeline_id, replace, rngseed):

	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["targyears"]
	scenario = my_data["scenario"]
	baseyear = my_data["baseyear"]
	ais_samples = my_data["ais_samps"]
	eais_samples = my_data["eais_samps"]
	wais_samples = my_data["wais_samps"]
	gis_samples = my_data["gis_samps"]

	# Generate the sample indices
	rng = np.random.default_rng(rngseed)
	sample_inds = rng.choice(ais_samples.shape[0], size=nsamps, replace=replace)

	# Store the samples for AIS components
	eais_samps = eais_samples[sample_inds,:]
	wais_samps = wais_samples[sample_inds,:]
	ais_samps = ais_samples[sample_inds,:]
	gis_samps = gis_samples[sample_inds,:]

	WriteOutput(eais_samps, wais_samps, ais_samps, gis_samps, years, scenario, baseyear, pipeline_id, nsamps)
	return(0)



def bamber19_project_icesheets_temperaturedriven(climate_data_file, pipeline_id, replace, rngseed):
	# Set rng
	rng = np.random.default_rng(rngseed)

	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["targyears"]
	scenario = my_data["scenario"]
	baseyear = my_data["baseyear"]
	ais_samplesH = my_data["ais_sampsH"]
	eais_samplesH = my_data["eais_sampsH"]
	wais_samplesH = my_data["wais_sampsH"]
	gis_samplesH = my_data["gis_sampsH"]
	ais_samplesL = my_data["ais_sampsL"]
	eais_samplesL = my_data["eais_sampsL"]
	wais_samplesL = my_data["wais_sampsL"]
	gis_samplesL = my_data["gis_sampsL"]

	# identify which samples to draw from high vs low scenarios
	useHigh=pickScenario(climate_data_file, scenario,rng);
	nsamps=useHigh.size


	# Generate the sample indices

	sample_inds = rng.choice(ais_samplesL.shape[0], size=nsamps, replace=replace)

	# Store the samples for AIS components
	eais_samps = eais_samplesL[sample_inds,:]
	wais_samps = wais_samplesL[sample_inds,:]
	ais_samps = ais_samplesL[sample_inds,:]
	gis_samps = gis_samplesL[sample_inds,:]

	eais_samps[useHigh,:] = eais_samplesH[sample_inds[useHigh],:]
	wais_samps[useHigh,:] = wais_samplesH[sample_inds[useHigh],:]
	ais_samps[useHigh,:] = ais_samplesH[sample_inds[useHigh],:]
	gis_samps[useHigh,:] = gis_samplesH[sample_inds[useHigh],:]

	WriteOutput(eais_samps, wais_samps, ais_samps, gis_samps, years, scenario, baseyear, pipeline_id, nsamps)
	return(0)

def WriteOutput(eais_samps, wais_samps, ais_samps, gis_samps, years, scenario, baseyear, pipeline_id, nsamps):
    # Store the variables in a pickle
	output = {'eais_samps': eais_samps, 'wais_samps': wais_samps, \
				'ais_samps': ais_samps, 'gis_samps': gis_samps, 'years': years, \
				'scenario': 'temperature-driven', 'baseyear': baseyear}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the projections to the netCDF files
	WriteNetCDF(pipeline_id, eais_samps, years, nsamps, "EAIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, wais_samps, years, nsamps, "WAIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, ais_samps, years, nsamps, "AIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, gis_samps, years, nsamps, "GIS", scenario, baseyear)

def WriteNetCDF(pipeline_id, global_samps, years, nsamps, ice_source, scenario, baseyear):

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_{1}_globalsl.nc".format(pipeline_id, ice_source))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(years))
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
	rootgrp.description = "Global SLR contribution from {0} from the Bamber et al. 2019 IPCC AR6 workflow".format(ice_source)
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}".format(pipeline_id)
	rootgrp.scenario = scenario
	rootgrp.baseyear = baseyear
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = years
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = global_samps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()


def GetSATData(fname, scenario, refyear_start=1850, refyear_end=1900, year_start=1900, year_end=2300):

	# Open the file
	df_ssp = h5py.File(fname,'r')

	# Extract the surface temperature for this scenario
	if scenario not in df_ssp.keys():
		raise Exception("Scenario {} not found in {}".format(scenario, fname))
	ssp_folder = df_ssp.get(scenario)
	sat_ssp = ssp_folder.get('surface_temperature')

	# Get the number of ensemble members from the data
	_,nens = sat_ssp.shape

	# Extract the years available
	sat_years = df_ssp['year'][()]

	# Which indices align with the reference and trim years
	refyear_start_idx = np.flatnonzero(sat_years == refyear_start)[0]
	refyear_end_idx = np.flatnonzero(sat_years == refyear_end)[0] #+ 1
	year_start_idx = np.flatnonzero(sat_years == year_start)[0]
	year_end_idx = np.flatnonzero(sat_years == year_end)[0] + 1

	# Normalize and crop temperature series
	Time = np.arange(year_start, year_end+1)
	SATave = np.mean(sat_ssp[refyear_start_idx:refyear_end_idx,:], axis=0)
	SAT = sat_ssp[year_start_idx:year_end_idx,:] - SATave

	# Close the h5 file
	df_ssp.close()

	# Done
	return(SAT, Time, nens)

def pickScenario(climate_data_file, scenario,rng):
	# Load the temperature data
	
	SAT,Time,NumTensemble = GetSATData(climate_data_file, scenario)

	# find integrated SAT over 2000-2099
	x2=np.where((Time[:]<2100) * (Time[:]>=2000))
	SAT2=SAT[x2]
	iSAT=SAT2.sum(axis=0)

	# convert integrated temperature into a normalized variable between low and high scenarios

	# Bamber 19 low scenario: 0.7 C in 2000, 1.5 C in 2050,
	# 2.0 C in 2100 = 70 + 20 + 40 + 12.5 = 142.5 C*yr if I've done this correctly

	# Bamber 19 high scenario: 0.7 C in 2000, 2.0 C in 2050, 5.0 C in 2100 =
	# 70 + 32.5 + 65 + 75 = 242.5 C*yr if I've done this correctly

	iSAT_marker = np.array([142.5,242.5])

	f2=np.minimum(1,np.maximum(0,(iSAT-iSAT_marker[0])/(iSAT_marker[1]-iSAT_marker[0])))
	weights=f2

	# Select which scenario to draw from for each sample
	selector = rng.random(iSAT.size)
	useHigh = (selector<weights)
	return useHigh

if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the IPCC AR6 Bamber et al. 2019 ice sheet projection stage.",\
	epilog="Note: This is meant to be run as part of the ipccar6 module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 10)", default=10, type=int)
	parser.add_argument('--replace', help="Allow sampling with replacement (default = 1)", choices=(0,1), type=int, default=1)
	parser.add_argument('--seed', help="Seed for the random number generator (default = 1234)", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str, default="")


	# Parse the arguments
	args = parser.parse_args()

	if len(args.climate_data_file) == 0:
		# Run the preprocessing stage with the user defined RCP scenario
		bamber19_project_icesheets(args.nsamps, args.pipeline_id, args.replace, args.seed)
	else:
		bamber19_project_icesheets_temperaturedriven(args.climate_data_file,args.pipeline_id, args.replace, args.seed)

	exit()
