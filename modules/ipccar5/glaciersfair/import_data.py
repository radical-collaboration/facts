import os
import sys
import re
import numpy as np

'''
import_data.py

Imports the emulated projections provided by Tamsin.

Parameters:
filename		Name of the file to import

Return:
data_dict		Dictionary of data using the header line as keys. The key for the data 
				is 'data'.

'''


def import_data(filename):
	
	# Initialize the data dictionary structure
	data_dict = {}
	
	# Load the emulated data file
	with open(filename, 'r') as f:
		
		# Get the header info
		header_line = f.readline().rstrip()
		header = np.array(header_line.split(","))
		
		# Up to what index is meta data?
		first_year_idx = np.flatnonzero(header == "1850")
		
		# Define the data indices
		data_idx = np.arange(first_year_idx,len(header))
		
		# Define the meta data indices
		meta_idx = np.arange(first_year_idx)
		
		# Extract the years from the header line
		data_dict['years'] = np.array([int(x) for x in header[data_idx]])
		
		# Initialize the data dictionary
		for i in meta_idx:
			data_dict[header[i]] = []
		data_dict['data'] = []
		
		# Loop through the lines in the file
		for line in f:
			
			# Get the line pieces
			line = line.rstrip()
			line_pieces = np.array(line.split(","))
			
			# Replace NA with numpy nan
			line_pieces[np.where(line_pieces == "NA")] = np.nan

			# Extract the data from the line
			data_dict['data'].append([float(x) for x in line_pieces[data_idx]])
			
			# Extract the meta data
			for i in meta_idx:
				data_dict[header[i]].append(line_pieces[i])
	
	# Convert everything into numpy arrays
	for this_key in data_dict.keys():
		data_dict[this_key] = np.array(data_dict[this_key])
	
	# Return the data dictionary
	return(data_dict)

	

if __name__ == "__main__":
	
	filename = "CLIMATE_FORCING_1850.csv"
	
	x = import_data(filename)
	
	print(x['data'].shape)
	
	exit()