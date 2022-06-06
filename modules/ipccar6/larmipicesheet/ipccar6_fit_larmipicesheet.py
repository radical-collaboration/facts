import os
import argparse
import pickle
import numpy as np

''' ipccar6_fit_larmipicesheet.py

This runs the fitting stage for the ISMIP Raw ice sheet component of the IPCC AR6
workflow.  Since there's no fitting done when performing a direct sample from a discrete
distribution, the fitting performed here is the fitting done for the surface mass balance
portion.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

Notes:
GMSLR projection program used for IPCC WG1 AR5
Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19
Adapted for use in FACTS by Gregory Garner 20 November 2019

'''

def ipccar6_fit_larmipicesheet(pipeline_id):

	# Read in the preprocessed data
	data_file = "{}_smbdata.pkl".format(pipeline_id)
	try:
		f = open(data_file, 'rb')
	except:
		print("Cannot open data file: {}\n".format(data_file))

	# Extract the data variables
	my_config = pickle.load(f)
	f.close()

	scenario = my_config["scenario"]

	# m SLE from Antarctica during 1996 to 2005 according to AR5 chapter 4
	dant = (2.37+0.13)*1e-3

	# m SLE from Greenland during 1996 to 2005 according to AR5 chapter 4
	dgreen = (3.21-0.30)*1e-3

	# Fraction of SLE from Greenland during 1996 to 2005 assumed to result from
	# rapid dynamical change, with the remainder assumed to result from SMB change
	fgreendyn = 0.5

	# Conversion factor for Gt to m SLE
	mSLEoGt = 1e12/3.61e14*1e-3

	# Greenland SMB ---------------------------------------------------------
	dtgreen=-0.146 # Delta_T of Greenland ref period wrt AR5 ref period
	fnlogsd=0.4 # random methodological error of the log factor
	febound=[1,1.15] # bounds of uniform pdf of SMB elevation feedback factor

	# Antarctic SMB ---------------------------------------------------------
	# The following are [mean,SD]
	pcoK=[5.1,1.5] # % change in Ant SMB per K of warming from G&H06
	KoKg=[1.1,0.2] # ratio of Antarctic warming to global warming from G&H06
	smax=0.35 # max value of S in 13.SM.1.5

	# Greenland Dynamic -----------------------------------------------------
	if scenario in ['rcp85','ssp5_85','ssp585']:
		gdyn_finalrange=[0.020,0.085]
	else:
		gdyn_finalrange=[0.014,0.063]

	# Antarctic Dynamic -----------------------------------------------------
	# For SMB+dyn during 2005-2010 Table 4.6 gives 0.41+-0.24 mm yr-1 (5-95% range)
	# For dyn at 2100 Chapter 13 gives [-20,185] mm for all scenarios
	adyn_finalrange = [-0.020,0.185]
	adyn_startratemean = 0.41
	adyn_startratepm = 0.20

	# Write the fitted parameters to a pickle file
	output = {'dant': dant, 'dgreen': dgreen, 'fgreendyn': fgreendyn, 'mSLEoGt': mSLEoGt, \
				'dtgreen': dtgreen, 'fnlogsd': fnlogsd, 'febound': febound, 'pcoK': pcoK, \
				'KoKg': KoKg, 'smax': smax, 'gdyn_finalrange': gdyn_finalrange, \
				'adyn_finalrange': adyn_finalrange, 'adyn_startratemean': adyn_startratemean, \
				'adyn_startratepm': adyn_startratepm}

	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(0)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the IPCC AR6 LARMIP ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the ipccar6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Run the preprocessing stage with the user defined RCP scenario
	ipccar6_fit_larmipicesheet(args.pipeline_id)

	# Done
	exit()
