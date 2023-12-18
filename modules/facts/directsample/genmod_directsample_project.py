import numpy as np
import argparse
import pickle
import os
import re

''' genmod_project_directsample.py

The "genmod - directsample" module ingests a 2-dimensional array that contains
samples of a sea-level change contributor (i.e. ice sheets, glaciers, etc.) as a function
of time.  These samples are treated as if they represent the true distributiuon of this
contributor over time.

This stage, as the name of the module implies, directly samples from the data provided in
the preprocessing stage.

The method of sampling depends on the number of samples requested. By default, the user 
may request up to the number of samples provided by the original data set. Asking for more 
would result in an error. 

Setting 'bootstrap' to 1 allows for replacement during the sampling process. This also
allows the user to request more samples than the original data set provides.

Parameters: 
nsamps              Number of samples to produce
replace             Allow sampling with replacement
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline


Output:
"{pipeline_id}_projections.pkl" = User requested samples

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def genmod_project_directsample(nsamps, pipeline_id, replace, rngseed):
	
	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)
	
	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["years"]
	samples = my_data["samples"]
	
	# Extract the pipeline id from the data file name
	#regex_match = re.search("genmod_directsample_(.*)_data\.pkl", datafile)
	#pipeline_id = regex_match.group(1)
	
	# Generate the sample indices
	rng = np.random.default_rng(rngseed)
	sample_inds = rng.choice(samples.shape[0], size=nsamps, replace=replace)
	
	# Store the samples
	samps = samples[sample_inds,:]
	   
    # Store the variables in a pickle
	output = {'samps': samps, 'years': years}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the generic module \'directsample\' projection stage.",\
	epilog="Note: This is meant to be run as part of the \'genmod\' module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 10)", default=10)
	parser.add_argument('--replace', help="Allow sampling with replacement (default = 0)", choices=('0','1'), default=0)
	parser.add_argument('--seed', help="Seed for the random number generator (default = 1234)", default=1234)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection stage with the provided arguments
	genmod_project_directsample(int(args.nsamps), args.pipeline_id, args.replace, int(args.seed))
	
	exit()