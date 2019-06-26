import numpy as np
import pickle
import sys
import os
import argparse
from scipy.stats import norm
from scipy.stats import t

''' kopp14_project_thermalexp.py

This script runs the thermal expansion projection task for the Kopp 2014 workflow. 
This task generates global contributions to sea-level change due to thermosteric heating
of the ocean.

Parameters: 
configfile = Pickle file containing the thermal expansion component's configuration.
			 This is produced by kopp14_preprocess_thermalexp.py.
fitfile = Pickle file produced by kopp14_fit_thermalexp.py
nsamps = Numer of samples to produce
seed = Seed for the random number generator

Output: Pickle file containing the global sea-level rise projections

'''

def kopp14_project_thermalexp(configfile, fitfile, nsamps, seed):
	
	# Load the configuration file
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	#rcp_scenario = my_config["rcp_scenario"]
	#datayears = my_config["datayears"]
	targyears = my_config["targyears"]
	#mergeZOSZOSTOGA = my_config["mergeZOSZOSTOGA"]
	#smoothwin = my_config["smoothwin"]
	#driftcorr = my_config["driftcorr"]
	#baseyear = my_config["baseyear"]
	GCMprobscale = my_config["GCMprobscale"]
	
	# Load the fit file
	try:
		f = open(fitfile, 'rb')
	except:
		print("Cannot open fit file\n")
	
	# Extract the fit variables
	my_fit = pickle.load(f)
	f.close()
	
	ThermExpMean = my_fit["ThermExpMean"]
	ThermExpStd = my_fit["ThermExpStd"]
	ThermExpYears = my_fit["ThermExpYears"]	
	ThermExpN = my_fit["ThermExpN"]
	
	# Evenly sample an inverse normal distribution and permutate it
	# Note: This may be a bug being ported over from Kopp 2014 which could result in 
	# 		overconfident projections
	np.random.seed(seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	norm_inv_perm = np.random.permutation(norm_inv)
	
	# Determine the scale coefficient
	ThermExpScale = norm.ppf(0.95)/norm.ppf(GCMprobscale)
	
	## Generate the samples --------------------------------------------------------------
	# Initialize variable to hold the samples
	thermsamps = np.full((nsamps, len(targyears)), np.nan)
	
	# Loop over the target years
	for i in np.arange(0,len(targyears)):
		
		# Find the index of ThermExp* that matches this target year
		this_year_ind = np.flatnonzero(ThermExpYears == targyears[i])
		
		# Generate the samples for this year
		#temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpN[this_year_ind]-1)
		temp = t.ppf(norm.cdf(norm_inv_perm), ThermExpN[i]-1)  # Replicates bug in K14 master code
		thermsamps[:,i] = (ThermExpScale * temp * ThermExpStd[this_year_ind]) + ThermExpMean[this_year_ind]

	# Save the global thermal expansion projections to a pickle
	output = {"thermsamps": thermsamps}
	outfile = open(os.path.join(os.path.dirname(__file__), "kopp14_thermexp_projections.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion projection stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', '-n', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', '-s', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	
	parser.add_argument('--config_file', help="Configuration file produced in the pre-processing stage",\
	default=os.path.join(os.path.dirname(__file__), "kopp14_thermalexp_config.pkl"))
	
	parser.add_argument('--fit_file', help="Fit file produced in the fitting stage",\
	default=os.path.join(os.path.dirname(__file__), "kopp14_thermalexp_fit.pkl"))
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	kopp14_project_thermalexp(args.config_file, args.fit_file, args.nsamps, args.seed)
	
	exit()