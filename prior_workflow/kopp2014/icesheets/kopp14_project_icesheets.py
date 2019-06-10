import pickle
import sys
import os
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
parfile = Name of the file containing the fitted distribution data
corfile = Name of the file containing the correlation information
seed = The seed to use for the random number gerneration

Output: Pickled file containing "nsamps" of global sea-level rise due to ice sheets

Note: The 'infile' must contain the variables 'batheteais', 'bathetwais', 'bathetgis', 
'arthetais', 'arthetgis', and 'islastdecade' within a single dictionary.

'''

def kopp14_project_icesheets(nsamps, parfile, corfile, seed=1234):
	
	## Read in the fitted parameters from parfile
	# Open the file
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
	outdir = os.path.join(os.path.dirname(__file__), "data")
	outfile = open(os.path.join(outdir, "kopp14_icesheets_projections.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	

if __name__ == '__main__':
	
	# Run the projection process for the parameters specified from the command line argument
	try:
		nsamps = int(sys.argv[1])
	except:
		nsamps = 20000
	
	try:
		fit_file = sys.argv[2]
	except:
		fit_file = os.path.join(os.path.dirname(__file__), "data", "kopp14_icesheets_fit.pkl")
	
	try:
		correlation_file = sys.argv[3]
	except:
		correlation_file = os.path.join(os.path.dirname(__file__), "data", "kopp14_icesheets_corr.pkl")
		
	try:
		seed = int(sys.argv[4])
	except:
		seed = 1234
	finally:
		kopp14_project_icesheets(nsamps, fit_file, correlation_file, seed)
	
	# Done
	exit()