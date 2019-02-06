import pickle
import sys
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
	outfile = open("kopp14_icesheets_fit.pkl", 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	

if __name__ == '__main__':
	
	# Run the fitting process on the file specified from the command line argument
	preprocess_file = sys.argv[1]
	kopp14_fit_icesheets(preprocess_file)
	
	# Done
	exit()