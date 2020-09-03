import numpy as np
import re

''' read_locationfile.py

Reads in the location file in order to get site names, site ids, lats, and lons

Parameters:
location_file = Background Rate file

'''

def ReadLocationFile(location_file):

	# Initialize variables to hold data and site information
	names = []
	ids = []
	lats = []
	lons = []
	
	# Compile the regex for finding commented lines
	comment_regex = re.compile(r'^#')
	
	# Open the rate file
	with open(location_file, 'r') as f:
		
		# Loop over the lines of the file
		for line in f:
			
			# Skip commented lines
			if re.search(comment_regex, line):
				continue
			
			# Split the line into components
			(this_name, this_id, this_lat, this_lon) = line.split("\t")

			# Store the information
			names.append(this_name)
			ids.append(int(this_id))
			lats.append(float(this_lat))
			lons.append(float(this_lon))
	
	# Cast everything as numpy arrays
	names = np.array(names)
	ids = np.array(ids)
	lats = np.array(lats)
	lons = np.array(lons)
	
	# Return variables
	return(names, ids, lats, lons)


if __name__ == "__main__":
	
	import argparse
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Reads in a location file used for localization in FACTS",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('infile', help="File from which to read input data")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Load the locations from the file
	(names, ids, lats, lons) = ReadLocationFile(args.infile)
	
	# Print the info to the screen
	print("--- infile: {} ---".format(args.infile))
	for i in np.arange(len(ids)):
		print("Name: {}\tID: {}\tLat: {}\tLon: {}".format(names[i], ids[i], lats[i], lons[i]))
	
	# Done
	exit()
	