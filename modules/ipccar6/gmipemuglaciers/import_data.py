import os
import sys
import re
import numpy as np


'''
import_data.py

Imports the emulated projections provided by Tamsin.

Parameters:
filename		Name of the file to import
model_driver	Name of the model set that drove the projections (CMIP6 or FAIR)

Return:
data_dict		Dictionary of data using the header line as keys. The key for the data 
				is 'data'.

'''


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
	
	filename = "191208_annual_CMIP6/projections_SSP585.csv"
	
	x = import_data(filename)
	
	print(x.keys())
	
	exit()