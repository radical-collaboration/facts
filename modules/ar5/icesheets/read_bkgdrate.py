import numpy as np
import pickle as p
import os
import sys
import argparse

''' read_bkgdrate.py

Reads in the background rate file in order to get site ids, lats, and lons

Parameters:
rate_file = Background Rate file
meta = Return just the site information (default: True)

'''

def read_bkgdrate(rate_file, meta=True):

	# Initialize variables to hold data and site information
	names = []
	ids = []
	lats = []
	lons = []
	rates = []
	sds = []
	
	# Define the rate file name
	#rate_file = os.path.join(os.path.dirname(__file__), "bkgdrate.tsv")
	
	# Open the rate file
	with open(rate_file, 'r') as f:
		
		# Skip the first line (header info)
		line = f.readline()
		
		# Loop over the lines of the file
		for line in f:
			
			# Split the line into components
			(this_name, this_id, this_lat, this_lon, this_rate, this_sd) = line.split("\t")
			
			# Store the information
			names.append(this_name)
			ids.append(int(this_id))
			lats.append(float(this_lat))
			lons.append(float(this_lon))
			rates.append(float(this_rate))
			sds.append(float(this_sd))
	
	# Cast everything as numpy arrays
	names = np.array(names)
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	rates = np.array(rates)
	sds = np.array(sds)
	
	# Return variables
	if(meta):
		return(names, ids, lats, lons)
	else:
		return(names, ids, lats, lons, rates, sds)
