import os
import sys
import re
import numpy as np

'''
filter_data

Filters Tamsin's emulated projections based on metadata fields

Parameters:
data_dict		Data dictionary that contains all the metadata and normal data
METADATA		Parameters named after the metadata fields (i.e. 'scenario', 'years', 
				'model', etc) that contain numpy arrays of said metadata field that you
				want to keep in the output

Return:
filtered_data_dict	Data dictionary filtered to user's request

'''

def filter_data(data_dict, model_driver, emulator_n=None, param_n=None, region=None, \
				GCM=None, year=None, ice_source=None, sample=None, scenario_sample=None):
	
	# Initialize the filtered data dictionary
	filtered_data_dict = {}
	
	# Initialize the index variables
	keep_idx = True
	
	# Filter the data------------------------------------------------
	# This section tests each input parameter to see if anything is there. If so, the
	# appropriate metadata is tested to see if it contains any of the items in the list.
	# If so, keep_idx are the indices to keep based on the filtered results.
	
	if year is not None:
		keep_idx *= np.isin(data_dict['year'], year)
	
	if region is not None:
		keep_idx *= np.isin(data_dict['region'], region)
	
	if GCM is not None:
		keep_idx *= np.isin(data_dict['GCM'], GCM)
	
	if ice_source is not None:
		keep_idx *= np.isin(data_dict['ice_source'], ice_source)
	
	# If this is CMIP6 data
	if(model_driver == "CMIP6"):
	
		if emulator_n is not None:
			keep_idx *= np.isin(data_dict['emulator_n'], emulator_n)
	
		if param_n is not None:
			keep_idx *= np.isin(data_dict['param_n'], param_n)
	
	# Otherwise, this is FAIR data
	else:
		
		if scenario_sample is not None:
			keep_idx *= np.isin(data_dict['scenario-sample'], scenario_sample)
		
		if sample is not None:
			keep_idx *= np.isin(data_dict['sample'], sample)
	
	
	# Get the index positions of the rows to keep
	keep_idx = np.flatnonzero(keep_idx)
	
	# If there are no matches, produce an error
	if not list(keep_idx):
		raise Exception("No matches for provided filter parameters")
	
	# If the user didn't pass anything to filter in the first place, return the original
	# data dictionary, but let the user know.
	if(len(keep_idx) <= 1):
		print("Nothing to filter. Returning original data dictionary")
		return(data_dict)
	
	# Loop through the keys of the dictionary
	for this_key in data_dict.keys():
			
		# Filter the metadata
		filtered_data_dict[this_key] = data_dict[this_key][keep_idx]
	
	# Done, return the filtered dictionary
	return(filtered_data_dict)


def import_data(filename, model_driver):
	
	# Initialize the data dictionary structure
	data_dict = {}
	
	# Load the emulated data file
	with open(filename, 'r') as f:
		
		# Get the header info
		header_line = f.readline().rstrip()
		header = np.array(header_line.split(","))
		meta_idx = np.arange(len(header))
		
		# Initialize the data dictionary
		for i in meta_idx:
			data_dict[header[i]] = []
		
		# Loop through the lines in the file
		for line in f:
			
			# Get the line pieces
			line = line.rstrip()
			line_pieces = np.array(line.split(","))
			
			# Replace NA with numpy nan
			line_pieces[np.where(line_pieces == "NA")] = np.nan
			
			# Extract the meta data
			for i in meta_idx:
				data_dict[header[i]].append(line_pieces[i])
	
	# Convert the data from string into float and from cm to mm
	data_dict['SLE'] = [float(x) * 10.0 for x in data_dict['SLE']]
	
	# If this is CMIP6 data
	if(model_driver == "CMIP6"):
	
		# Convert the year, param_n, and emulator_n from string into int
		data_dict['year'] = [int(x) for x in data_dict['year']]
		data_dict['param_n'] = [int(x) for x in data_dict['param_n']]
		data_dict['emulator_n'] = [int(x) for x in data_dict['emulator_n']]
	
	# Otherwise this is FAIR data
	else:
	
		# Convert the year and sample from string into int
		data_dict['year'] = [int(x) for x in data_dict['year']]
		data_dict['sample'] = [int(x) for x in data_dict['sample']]
		
	
	# Convert everything into numpy arrays
	for this_key in data_dict.keys():
		data_dict[this_key] = np.array(data_dict[this_key])
	
	# Return the data dictionary
	return(data_dict)


if __name__ == "__main__":
	
	from netCDF4 import Dataset
	import time
	
	# File name and driver
	scenario = "SSP585"
	filename = "191220_emulated/projections_FAIR_{}.csv".format(scenario)
	model_driver = "FAIR"
	
	# Import data
	x = import_data(filename, model_driver)
	
	# Initialize data structure
	glac_years = np.arange(2020,2101,10)
	glac_samples = np.arange(500) + 1
	glacier_data = np.zeros((len(glac_samples), len(glac_years)))
	
	# Loop through the filter elements
	for i in np.arange(len(glac_years)):
		
		print("Processing samples for year {}".format(glac_years[i]))
		
		for j in np.arange(len(glac_samples)):
			
			# Filter the data
			x_filtered = filter_data(x, model_driver=model_driver, ice_source="Glaciers", year=glac_years[i], sample=glac_samples[j])
			
			# Put this data into the data structure
			glacier_data[j,i] = np.sum(x_filtered['SLE'])
	
	# Write to netcdf file
	rootgrp = Dataset(os.path.join(os.path.dirname(__file__), "emulated_glaciers_{}.nc".format(scenario)), "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(glac_years))
	samp_dim = rootgrp.createDimension("samples", len(glac_samples))

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i4", ("samples",))

	# Create a data variable
	glac_sle = rootgrp.createVariable("glacier_sle", "f4", ("samples", "years"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Glaciers and ice caps from Tamsin's emulator"
	rootgrp.history = "Created " + time.ctime(time.time())
	glac_sle.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = glac_years
	samp_var[:] = glac_samples
	glac_sle[:,:] = glacier_data

	# Close the netcdf
	rootgrp.close()
	
	exit()
