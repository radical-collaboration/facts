import pickle
import sys
import os
import argparse
import numpy as np
from SampleISDists import SampleISDists
from ProjectGSL import ProjectGSL
import cholcov

''' kopp14_project_icesheets.py

This script runs the ice sheet projection task for the Kopp 2014 workflow. This task
requires data produced from the pre-processing and fit/calibration stages. This task
also requires at the time of call the number of samples to produce and the names of the
input files generated from the previous stages.

Parameters: 
nsamps = Number of samples to produce
seed = The seed to use for the random number gerneration
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickled file containing "nsamps" of global sea-level rise due to ice sheets

'''

def kopp14_project_icesheets(nsamps, seed, pipeline_id):
	
	## Read in the fitted parameters from parfile
	# Open the file
	parfile = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(parfile, 'rb')
	except:
		print("Cannot open parfile\n")
	
	# Extract the data from the file
	my_data = pickle.load(f)
	batheteais = my_data['batheteais']
	bathetwais = my_data['bathetwais']
	bathetgis = my_data['bathetgis']
	arthetais = my_data['arthetais']
	arthetgis = my_data['arthetgis']
	islastdecade = my_data['islastdecade']
	f.close()
	
	## Read in the correlation information
	corfile = "{}_corr.pkl".format(pipeline_id)
	try:
		f = open(corfile, 'rb')
	except:
		print("Cannot open corfile\n")
	
	# Extract the data
	my_data = pickle.load(f)
	bacorris = my_data['bacorris']
	arcorris = my_data['arcorris']
	f.close()
	
	# Generate samples of ice sheet accelerations
	sigmas = np.array([bathetgis[2], bathetwais[2], batheteais[2]])
	mus = np.array([bathetgis[1], bathetwais[1], batheteais[1]])
	offsets = np.array([bathetgis[0], bathetwais[0], batheteais[0]])
	baissamps = SampleISDists(nsamps, sigmas, mus, offsets, islastdecade, bacorris, seed)
	
	sigmas = np.array([arthetgis[2], arthetais[2]])
	mus = np.array([arthetgis[1], arthetais[1]])
	offsets = np.array([arthetgis[0], arthetais[0]])
	arislastdecade = np.array([islastdecade[0], islastdecade[1]+islastdecade[2]])
	arissamps = SampleISDists(nsamps, sigmas, mus, offsets, arislastdecade, arcorris, seed+1234)
	
	# Project global sea-level rise over time
	targyears = np.arange(2010, 2101, 10)
	(arsamps, basamps, hysamps) = ProjectGSL(baissamps, arissamps, islastdecade, targyears)
	
	# Put the results into a dictionary
	output = {'arsamps': arsamps, 'basamps': basamps, 'hysamps': hysamps}
	
	# Write the results to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process for the parameters specified from the command line argument
	kopp14_project_icesheets(args.nsamps, args.seed, args.pipeline_id)
	
	# Done
	exit()