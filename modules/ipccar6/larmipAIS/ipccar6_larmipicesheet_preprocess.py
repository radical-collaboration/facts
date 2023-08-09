import argparse
import os
import sys
import re
import pickle as p
import numpy as np
import fnmatch
from netCDF4 import Dataset
from Import2lmData import *

''' ipccar6_preprocess_larmipicesheet.py

This runs the preprocessing stage for the LARMIP ice sheet component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code
scenario = RCP scenario of interest

'''

def ipccar6_preprocess_larmipicesheet(pipeline_id, pyear_start, pyear_end, pyear_step, baseyear, scenario):

	# Get the searchable RCP scenario format
	scenario_dict = {'rcp85': "RCP85", 'rcp70': "RCP70", 'rcp60': "RCP70", 'rcp45': "RCP45", 'rcp26': "RCP26", 'rcp19': "RCP19", \
						'ssp119': "RCP19", 'ssp126': "RCP26", 'ssp245': "RCP45", 'ssp370': "RCP70", \
						'ssp434': "RCP45", 'ssp460': "RCP70", 'ssp585': "RCP85"}
	try:
		scenario_search = scenario_dict[scenario]
	except:
		scenario_search = scenario

	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Set the base year
	#baseyear = 2005

	# Initialize the data structures
	eais_samps = []
	wais_samps = []
	pen_samps = []
	ais_samps = []

	# Get a list of files for this scenario
	datadir = os.path.join(os.path.dirname(__file__), "larmip_data")	# Directory for FACTS module
	#datadir = os.path.join(os.path.dirname(__file__),"data", "larmip_data_original")	# Directory for testing
	filenames = fnmatch.filter(os.listdir(datadir), "*{}.nc".format(scenario_search))

	# Loop through the input files
	for this_filename in filenames:

		# Open the netCDF file
		nc = Dataset(os.path.join(datadir, this_filename), 'r')

		# Get the years available
		nc_years = nc.variables['Time'][:]

		# Get the base year values
		# Note: Assumes WAIS is the sum of Amundsen, Ross, and Weddell regions
		# Note: Instead of loading each of the other variables, I subtract EAIS and PEN
		# from AIS to get WAIS. This matches the sum of the other regions from my tests.
		baseyear_idx = np.flatnonzero(nc_years == baseyear)
		baseyear_eais_val = nc.variables['EAIS'][1:,baseyear_idx]
		baseyear_ais_val = nc.variables['Antarctica'][1:,baseyear_idx]
		baseyear_pen_val = nc.variables['Peninsula'][1:,baseyear_idx]
		baseyear_wais_val = baseyear_ais_val - baseyear_eais_val - baseyear_pen_val

		# Linear extrapolation to year 2100
		eais_2100 = (nc.variables['EAIS'][1:,-1] - nc.variables['EAIS'][1:,-2]) + nc.variables['EAIS'][1:,-1]
		ais_2100 = (nc.variables['Antarctica'][1:,-1] - nc.variables['Antarctica'][1:,-2]) + nc.variables['Antarctica'][1:,-1]
		pen_2100 = (nc.variables['Peninsula'][1:,-1] - nc.variables['Peninsula'][1:,-2]) + nc.variables['Peninsula'][1:,-1]
		#wais_2100 = ais_2100 - eais_2100 - pen_2100

		# Get the contributions for each target year
		nc_year_idx = np.isin(nc_years, targyears)
		eais_samps_temp = nc.variables['EAIS'][1:,nc_year_idx]
		ais_samps_temp = nc.variables['Antarctica'][1:,nc_year_idx]
		pen_samps_temp = nc.variables['Peninsula'][1:,nc_year_idx]

		# Append the year 2100 values
		eais_samps_temp = np.concatenate((eais_samps_temp, eais_2100[:,np.newaxis]), axis=1)
		ais_samps_temp = np.concatenate((ais_samps_temp, ais_2100[:,np.newaxis]), axis=1)
		pen_samps_temp = np.concatenate((pen_samps_temp, pen_2100[:,np.newaxis]), axis=1)

		# Calculate the WAIS
		wais_samps_temp = ais_samps_temp - eais_samps_temp - pen_samps_temp

		# Center these samples on the base year
		eais_samps_temp = eais_samps_temp - baseyear_eais_val
		ais_samps_temp = ais_samps_temp - baseyear_ais_val
		pen_samps_temp = pen_samps_temp - baseyear_pen_val
		wais_samps_temp = wais_samps_temp - baseyear_wais_val

		# Append these to the larger data structure
		eais_samps.append(eais_samps_temp)
		ais_samps.append(ais_samps_temp)
		pen_samps.append(pen_samps_temp)
		wais_samps.append(wais_samps_temp)

		# Close the netCDF file
		nc.close()

	# Convert the data structures into numpy arrays and convert sle from m to mm
	eais_samps = np.array(eais_samps) * 1000
	ais_samps = np.array(ais_samps) * 1000
	pen_samps = np.array(pen_samps) * 1000
	wais_samps = np.array(wais_samps) * 1000

	# Populate the output dictionary
	outdata = {'eais_samps': eais_samps, 'wais_samps': wais_samps, 'pen_samps': pen_samps, \
				'ais_samps': ais_samps, 'scenario': scenario, 'targyears': targyears, 'baseyear': baseyear}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()

	return(0)


'''
ipccar6_preprocess_larmipicesheet2lm()

Uses the larmip data produced from the 2 layer model data.

'''
def ipccar6_preprocess_larmipicesheet2lm(pipeline_id, pyear_start, pyear_end, pyear_step, baseyear, scenario):

	# Get the searchable RCP scenario format
	scenario_dict = {'rcp85': "SSP585", 'rcp70': "SSP370", 'rcp60': "SSP370", 'rcp45': "SSP245", 'rcp26': "SSP126", 'rcp19': "SSP119", \
						'ssp119': "SSP119", 'ssp126': "SSP126", 'ssp245': "SSP245", 'ssp370': "SSP370", \
						'ssp434': "SSP245", 'ssp460': "SSP370", 'ssp585': "SSP585"}

	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Set the base year
	#baseyear = 2005

	# Define the temperature limit scenarios that are available
	tlim_scenarios = ["ssp119", "ssp126", "ssp245", "ssp370", "ssp585"]

	# Initialize the data structures
	make_data_structure = True

	# Get a list of files for this scenario
	datadir = os.path.join(os.path.dirname(__file__), "larmip_data_2lm")	# Directory for FACTS module
	#datadir = os.path.join(os.path.dirname(__file__),"data", "larmip_data_2lm")	# Directory for testing

	# Get the projection indices for this run
	(_, idx_dict) = Import2lmData(scenario=scenario, directory=os.path.join(datadir, "..", "smb_data"), tlim_scenarios=tlim_scenarios)

	# Loop through the matching scenarios
	for this_scenario in idx_dict.keys():

		# Initialize arrays to hold samples from all the models
		eais_model_samps = []
		wais_model_samps = []
		pen_model_samps = []
		ais_model_samps = []

		# If there are no matching indices for this scenario, skip to the next
		if idx_dict[this_scenario] is None:
			continue

		# Get the list of model files that match this scenario
		scenario_search = scenario_dict[this_scenario]
		filenames = fnmatch.filter(os.listdir(datadir), "*{}.nc".format(scenario_search))

		# Loop through the input files
		for this_filename in filenames:

			# Open the netCDF file
			nc = Dataset(os.path.join(datadir, this_filename), 'r')

			# Turn off auto-masking so that we can properly test for missing data
			nc.set_auto_mask(False)

			# Get the years available
			nc_years = nc.variables['Time'][:]

			# Get the TLM indices from this model file
			tlm_idx = nc.variables['tlm_idx'][:]

			# Read in the samples
			eais_vals = nc.variables['EAIS'][...].astype(float)
			ais_vals = nc.variables['Antarctica'][...].astype(float)
			pen_vals = nc.variables['Peninsula'][...].astype(float)

			# Set missing data to np.nan
			eais_vals[eais_vals < -999] = np.nan
			ais_vals[ais_vals < -999] = np.nan
			pen_vals[pen_vals < -999] = np.nan

			wais_vals = ais_vals - eais_vals - pen_vals

			# Get the base year values
			# Note: Assumes WAIS is the sum of Amundsen, Ross, and Weddell regions
			# Note: Instead of loading each of the other variables, I subtract EAIS and PEN
			# from AIS to get WAIS. This matches the sum of the other regions from my tests.
			baseyear_idx = np.flatnonzero(nc_years == baseyear)
			baseyear_eais_val = eais_vals[:,baseyear_idx]
			baseyear_ais_val = ais_vals[:,baseyear_idx]
			baseyear_pen_val = pen_vals[:,baseyear_idx]
			baseyear_wais_val = baseyear_ais_val - baseyear_eais_val - baseyear_pen_val

			# Center these samples on the base year
			eais_vals = eais_vals - baseyear_eais_val
			ais_vals = ais_vals - baseyear_ais_val
			pen_vals = pen_vals - baseyear_pen_val
			wais_vals = wais_vals - baseyear_wais_val

			# Determine which years and samples to extract
			nc_year_idx = np.flatnonzero(np.isin(nc_years, targyears))[None,:]
			sample_idx = np.flatnonzero(np.isin(tlm_idx, idx_dict[this_scenario][:]))[:,None]

			# Initialize temporary arrays
			# This makes sure that ragged arrays from tlim scenarios don't cause problems
			eais_samps_temp = np.full((20000,nc_year_idx.shape[1]), np.nan)
			ais_samps_temp = np.full((20000,nc_year_idx.shape[1]), np.nan)
			pen_samps_temp = np.full((20000,nc_year_idx.shape[1]), np.nan)

			# Extract the samples
			eais_samps_temp[np.arange(len(sample_idx)),:] = eais_vals[sample_idx,nc_year_idx]
			ais_samps_temp[np.arange(len(sample_idx)),:] = ais_vals[sample_idx,nc_year_idx]
			pen_samps_temp[np.arange(len(sample_idx)),:] = pen_vals[sample_idx,nc_year_idx]

			#eais_samps_temp = eais_vals[sample_idx,nc_year_idx]
			#ais_samps_temp = ais_vals[sample_idx,nc_year_idx]
			#pen_samps_temp = pen_vals[sample_idx,nc_year_idx]

			# Calculate the WAIS
			wais_samps_temp = ais_samps_temp - eais_samps_temp - pen_samps_temp

			# Append these to the larger data structure
			eais_model_samps.append(eais_samps_temp)
			ais_model_samps.append(ais_samps_temp)
			pen_model_samps.append(pen_samps_temp)
			wais_model_samps.append(wais_samps_temp)

			# Close the netCDF file
			nc.close()

		# Convert model samples to numpy array
		eais_model_samps = np.array(eais_model_samps)
		ais_model_samps = np.array(ais_model_samps)
		pen_model_samps = np.array(pen_model_samps)
		wais_model_samps = np.array(wais_model_samps)

		# Create the uber structure for the samples or append to it
		if make_data_structure:
			eais_samps = eais_model_samps
			ais_samps = ais_model_samps
			pen_samps = pen_model_samps
			wais_samps = wais_model_samps
			make_data_structure = False
		else:
			eais_samps = np.concatenate((eais_samps, eais_model_samps), axis=1)
			ais_samps = np.concatenate((ais_samps, ais_model_samps), axis=1)
			pen_samps = np.concatenate((pen_samps, pen_model_samps), axis=1)
			wais_samps = np.concatenate((wais_samps, wais_model_samps), axis=1)

	# Convert the data structures into numpy arrays [models, samples, years]
	eais_samps = np.array(eais_samps)
	ais_samps = np.array(ais_samps)
	pen_samps = np.array(pen_samps)
	wais_samps = np.array(wais_samps)

	# Populate the output dictionary
	outdata = {'eais_samps': eais_samps, 'wais_samps': wais_samps, 'pen_samps': pen_samps, \
				'ais_samps': ais_samps, 'scenario': scenario, 'targyears': targyears, 'baseyear': baseyear}

	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()

	return(0)




'''
ipccar6_preprocess_ar5aissmb()

Pre-processes the data for the Antarctic surface mass balance
portion using the AR5 data and methodology.

'''

def ipccar6_preprocess_ar5aissmb(scenario, startyr, pipeline_id, pyear_start, pyear_end, pyear_step, tlm_flag):

	# Get the searchable RCP scenario format
	scenario_dict = {'rcp85': "rcp85", 'rcp70': "ssp370", 'rcp60': "rcp60", 'rcp45': "rcp45", 'rcp26': "rcp26", 'rcp19': "ssp119", \
						'ssp119': "ssp119", 'ssp126': "ssp126", 'ssp245': "ssp245", 'ssp370': "ssp370", \
						'ssp434': "ssp434", 'ssp460': "ssp460", 'ssp585': "ssp585"}
	try:
		scenario_search = scenario_dict[scenario]
	except:
		scenario_search = scenario

	# Define the target projection years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

	# Define the input data directory
	indir = os.path.join(os.path.dirname(__file__), "smb_data")	# Directory for FACTS module
	#indir = os.path.join(os.path.dirname(__file__), "data", "smb_data")		# Directory for testing

	# Load the two-layer model data
	if tlm_flag:

		# Import the data
		tlm_dict, _ = Import2lmData("surface_temperature", scenario_search, indir, refyear_start=1986, refyear_end=2005)

		# Filter the data for the appropriate years
		filtered_data_dict = Filter2lmData(tlm_dict, filter_years=np.arange(startyr,pyear_end+1))

		# Extract the years
		data_years = filtered_data_dict['years']

		# FOR THE TEMPORARY 2LM DATA ONLY - USES ONLY 44 UNIQUE TRAJECTORIES
		# Find the unique temperature trajectories
		#temp_samples = np.unique(filtered_data_dict["samples"], axis=0)
		temp_samples = filtered_data_dict["samples"]

		# Find the mean and sd of the ensemble
		temp_mean = np.nanmean(temp_samples, axis=0)
		temp_sd = np.nanstd(temp_samples, axis=0)

	else:

		# Define the input data files
		temp_mean_filename = "{0}_temperature_mean.nc".format(scenario_search)
		temp_sd_filename = "{0}_temperature_sd.nc".format(scenario_search)

		# Extract the data ------------------------------------------
		nc = Dataset(os.path.join(indir,temp_mean_filename), 'r')
		temp_years = 2006 + nc.variables['bound_time'][:,1]/360
		year_idx = np.flatnonzero(temp_years >= startyr)
		temp_mean = nc.variables['air_temperature'][year_idx]
		nc.close()

		nc = Dataset(os.path.join(indir,temp_sd_filename), 'r')
		temp_sd = nc.variables['air_temperature'][year_idx]
		nc.close()

		# Set the output data years to the temperature years
		data_years = temp_years

		# No direct temperature samples available.  Set to None.
		temp_samples = None

	# Integrate temperature to obtain K yr at ends of calendar years
	# Note - The original code I believe performs a cumulative sum of the standard
	# deviations, which is not correct.  Below I provide a fix to that bug as well as
	# a replication of the bug for diagnostic purposes.
	# Note - JG makes an assumption here so that standard deviations are able to be summed
	# over time.
	inttemp_mean = np.cumsum(temp_mean)
	#inttemp_sd = np.sqrt(np.cumsum(temp_sd**2))  # Fix the bug
	inttemp_sd = np.cumsum(temp_sd)  # Replicate the bug

	# Store preprocessed data in pickles
	output = {'temp_mean': temp_mean, 'temp_sd': temp_sd, 'inttemp_mean': inttemp_mean, \
				'inttemp_sd': inttemp_sd, 'data_years': data_years, 'startyr': startyr, \
				'scenario': scenario, 'targyears': targyears, 'temp_samples': temp_samples}

	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_smbdata.pkl".format(pipeline_id)), 'wb')
	p.dump(output, outfile)
	outfile.close()

	return(0)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the IPCC AR6 LARMIP ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the IPCC AR6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Scenario choices
	#scenario_choices = ['rcp85', 'rcp70', 'rcp60', 'rcp45', 'rcp26', 'rcp19',\
	#					'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp434', 'ssp460', 'ssp585']

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--scenario', help="Scenario of interest", default='rcp85')
	parser.add_argument('--baseyear', help="Year from which to start integrating temperature [default=2006]", type=int, default=2006)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--tlm_data', help="Use the two-layer model data [default=1,  use 2lm data]", default=1, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	if args.tlm_data == 0:
		# If projection year ends after 2100, stop because there's no larmip data beyond 2100
		if args.pyear_end > 2100:
			raise Exception("Projection end year is beyond 2100. pyear_end = {}".format(pyear_end))
		ipccar6_preprocess_larmipicesheet(args.pipeline_id, args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario)
	else:
		# If projection year ends after 2100, stop because there's no larmip data beyond 2100
		if args.pyear_end > 2300:
			raise Exception("Projection end year is beyond 2300. pyear_end = {}".format(pyear_end))
		ipccar6_preprocess_larmipicesheet2lm(args.pipeline_id, args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario)

	ipccar6_preprocess_ar5aissmb(args.scenario, args.baseyear, args.pipeline_id, args.pyear_start, args.pyear_end, args.pyear_step, args.tlm_data)

	# Done
	exit()
