import pickle
import sys
import os
import numpy as np
from CalcISDists import CalcISDists as calcDists

''' kopp14_fit_icesheets.py

This script runs the ice sheet fitting task for the Kopp 2014 workflow. Data produced
from the pre-processing task 'kopp14_preprocess_icesheets.py' are used here. This task 
requires a single paramter, which is the name of the pickled file produced by the 
pre-processing task. Output is another pickled file that contains the fitted parameters
for the log-normal distributions representing ice sheet contributions to future sea-level 
rise

Parameters: 
infile = Name of the file containing the pre-processed data

Output: Pickled file containing the fitted parameter values

Note: The 'infile' must contain the variables 'barates', 'lastdecadegt', and 'aris2090' 
within a single dictionary.

'''

def kopp14_fit_icesheets(infile):
	
	# Open the file
	try:
		f = open(infile, 'rb')
	except:
		print("Cannot open infile\n")
	
	# Extract the data from the file
	my_data = pickle.load(f)
	f.close()
	
	# Fit the distributions
	(batheteais, bathetwais, bathetgis, arthetais, arthetgis, islastdecade) = calcDists(my_data['barates'], my_data['lastdecadegt'], my_data['aris2090'])
	
	# Put the results into a dictionary
	output = {'batheteais': batheteais, 'bathetwais': bathetwais, 'bathetgis': bathetgis,\
				'arthetais': arthetais, 'arthetgis': arthetgis,\
				'islastdecade': islastdecade}
	
	# Write the results to a file
	outdir = os.path.join(os.path.dirname(__file__), 'data')
	outfile = open(os.path.join(outdir, "kopp14_icesheets_fit.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	

if __name__ == '__main__':
	
	# Run the fitting process using the provided rate file from preprocessing stage
	try:
		rate_file = sys.argv[1]
	except:
		rate_file = os.path.join(os.path.dirname(__file__), "data", "kopp14_icesheets_rates.pkl")
	finally:
		kopp14_fit_icesheets(rate_file)
	
	# Done
	exit()