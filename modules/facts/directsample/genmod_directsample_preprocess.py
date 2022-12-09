import numpy as np
import csv
import argparse
import pickle
import os

''' genmod_preprocess_directsample.py

The "genmod - directsample" module ingests a 2-dimensional array that contains
samples of a sea-level change contributor (i.e. ice sheets, glaciers, etc.) as a function
of time.  These samples are treated as if they represent the true distributiuon of this
contributor over time.

Preprocessing ingests the data from a provided input file and parses the data into a
format that can be sampled in a later stage.  The input file should be an ascii file with 
each row containing a single time-series sample (thus columns are the time-dimension).
The first row (i.e. the header) should contain the year associated with the samples in
that column.

Parameters: 
input_data_file     Name of the input data file
pipeline_id         Unique identifier to attach to this pipeline


Output:
"{pipeline_id}_data.pkl" = Data to be passed through pipeline

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def CountRowCol(file):
	with open(file, 'rU') as fp:
		x = csv.reader(fp)
		ncol = len(x.next())
		nrow = sum(1 for row in x) + 1
	return((nrow, ncol))

def genmod_preprocess_directsample(filename, pipeline_id):
	
	# Count the number of lines in the input file
	(nrow, ncol) = CountRowCol(filename)
	
	# Extract the data
	samples = np.empty((nrow-1, ncol))
	
	with open(filename, 'rU') as fp:
		reader = csv.reader(fp)
		i = 0
		years = [int(x) for x in reader.next()]
		for row in reader:
			samples[i,:] = row
			i += 1
    
    ###################################################
    # Store the data in a pickle
	output = {'years': years, 'samples': samples}
	
	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfilename = "{}_data.pkl".format(pipeline_id)
	outfile = open(os.path.join(outdir, outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the generic module \'directsample\' pre-processing stage.",\
	epilog="Note: This is meant to be run as part of the \'genmod\' module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--input_data_file', help="Full path and name of file to import")
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	genmod_preprocess_directsample(args.input_data_file, args.pipeline_id)
	
	exit()