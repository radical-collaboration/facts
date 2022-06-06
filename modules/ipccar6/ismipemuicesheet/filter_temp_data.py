import os
import sys
import re
import numpy as np

'''
filter_data

Filters Tamsin's climate forcing based on metadata fields

Parameters:
data_dict		Data dictionary that contains all the metadata and normal data
METADATA		Parameters named after the metadata fields (i.e. 'scenario', 'years', 
				'model', etc) that contain numpy arrays of said metadata field that you
				want to keep in the output

Return:
filtered_data_dict	Data dictionary filtered to user's request

'''

def filter_temp_data(data_dict, years=None, ensemble=None, GCM=None, scenario=None):
	
	# Initialize the filtered data dictionary
	filtered_data_dict = {}
	
	# Temporary variables to hold year-related data
	temp_years = data_dict['years']
	temp_data = data_dict['data']
	
	# Initialize the index variables
	keep_idx = True
	year_idx = None
	
	# Filter the data------------------------------------------------
	# This section tests each input parameter to see if anything is there. If so, the
	# appropriate metadata is tested to see if it contains any of the items in the list.
	# If so, keep_idx are the indices to keep based on the filtered results.
	
	if years is not None:
		year_idx = np.flatnonzero(np.isin(data_dict['years'], years))
		temp_years = temp_years[year_idx]
		temp_data = temp_data[:,year_idx]
	
	if scenario is not None:
		keep_idx *= np.isin(data_dict['scenario'], scenario)
	
	if GCM is not None:
		keep_idx *= np.isin(data_dict['GCM'], GCM)
	
	if ensemble is not None:
		keep_idx *= np.isin(data_dict['ensemble'], ensemble)
	
	# Get the index positions of the rows to keep
	keep_idx = np.flatnonzero(keep_idx)
	
	# If there are no matches, produce an error
	if not list(keep_idx):
		raise Exception("No matches for provided filter parameters")
	
	# If the user didn't pass anything to filter in the first place, return the original
	# data dictionary, but let the user know.
	if(len(keep_idx) <= 1 and year_idx is None):
		print("Nothing to filter. Returning original data dictionary")
		return(data_dict)
	
	# Loop through the keys of the dictionary
	for this_key in data_dict.keys():
		
		# Skip the 'years' and 'data' since they are different shapes than the rest of
		# the dictionary entries and are handled separately
		if (this_key == "years" or this_key == "data"):
			continue
		
		# Filter the metadata
		filtered_data_dict[this_key] = data_dict[this_key][keep_idx]
	
	# Filter the 'years' and 'data' entries
	filtered_data_dict['years'] = temp_years
	filtered_data_dict['data'] = temp_data[keep_idx,:]
		
	# Done, return the filtered dictionary
	return(filtered_data_dict)



def import_temp_data(filename):
	
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
	
	x = import_temp_data(filename)
	
	x_filtered = filter_temp_data(x, ensemble="FAIR", scenario=["SSP119", "SSP126", "SSP245", "SSP370", "SSP585"], years=2100)
	
	print(x_filtered['data'].shape)
	
	exit()
