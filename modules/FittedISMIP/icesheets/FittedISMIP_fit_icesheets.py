import argparse
import os
import sys
import pickle
import numpy as np

''' FittedISMIP_fit_icesheets.py

Fitting process for the ismip6gis icesheet module.

Fitting is done external to this FACTS workflow.  Parameters from this fitting process
are read in and stored here for the next stage.

Note: 'pipeline_id' is a unique identifier that distinguishes it among other instances
of this module within the same workflow.

'''

def FittedISMIP_fit_icesheets(pipeline_id):
	
	# Initialize dictionaries to hold model parameters
	groups_dict = {}
	models_dict = {}
	betas_dict = {}
	sigmas_dict = {}
	
	# Load the Greenland ice sheet model fits
	parmfile = os.path.join(os.path.dirname(__file__), "FittedParms_GrIS_ALL.csv")
	groups, models, betas, sigmas = ReadParameterFile(parmfile)
	groups_dict["GIS"] = groups
	models_dict["GIS"] = models
	betas_dict["GIS"] = betas
	sigmas_dict["GIS"] = sigmas
	
	# Load the West Antarctic ice sheet model fits
	parmfile = os.path.join(os.path.dirname(__file__), "FittedParms_AIS_WAIS.csv")
	groups, models, betas, sigmas = ReadParameterFile(parmfile)
	groups_dict["WAIS"] = groups
	models_dict["WAIS"] = models
	betas_dict["WAIS"] = betas
	sigmas_dict["WAIS"] = sigmas
	
	# Load the East Antarctic ice sheet model fits
	parmfile = os.path.join(os.path.dirname(__file__), "FittedParms_AIS_EAIS.csv")
	groups, models, betas, sigmas = ReadParameterFile(parmfile)
	groups_dict["EAIS"] = groups
	models_dict["EAIS"] = models
	betas_dict["EAIS"] = betas
	sigmas_dict["EAIS"] = sigmas
	
	# Load the Antarctic Peninsula ice sheet model fits
	parmfile = os.path.join(os.path.dirname(__file__), "FittedParms_AIS_PEN.csv")
	groups, models, betas, sigmas = ReadParameterFile(parmfile)
	groups_dict["PEN"] = groups
	models_dict["PEN"] = models
	betas_dict["PEN"] = betas
	sigmas_dict["PEN"] = sigmas
	
	# Define the linear trend terms to be added to the samples
	#trend_mean = {"EAIS": -0.02, "WAIS": 0.28, "PEN": 0.06, "GIS": 0.46}	# SOD
	#trend_sd = {"EAIS": 0.05, "WAIS": 0.03, "PEN": 0.01, "GIS": 0.04}		# SOD
	trend_mean = {"EAIS": 0.09, "WAIS": 0.18, "PEN": 0.06, "GIS": 0.19}
	trend_sd = {"EAIS": 0.04, "WAIS": 0.09, "PEN": 0.03, "GIS": 0.1}

	 ###################################################
    # Store the data in a pickle
	output = {'groups_dict': groups_dict, 'models_dict': models_dict, \
				'betas_dict': betas_dict, 'sigmas_dict': sigmas_dict, \
				'trend_mean': trend_mean, 'trend_sd': trend_sd}
	
	# Write the data to a file
	outdir = os.path.dirname(__file__)
	outfilename = "{}_fit.pkl".format(pipeline_id)
	outfile = open(os.path.join(outdir, outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	return(None)


def ReadParameterFile(file):
	
	# Initialize the return variables
	groups = []
	models = []
	betas = []
	sigmas = []
	
	# Open the file for reading
	with open(file, "r") as f:
		
		# Loop through the lines
		for line in f.readlines():
			
			# Strip off the endline character
			line = line.rstrip()
			
			# Split the line by commas
			line_pieces = line.split(",")
			
			# Put these data into the return arrays
			groups.append(line_pieces[0])
			models.append(line_pieces[1])
			sigmas.append(float(line_pieces[-1]))
			betas.append([float(x) for x in line_pieces[2:-1]])
	
	# Convert return arrays into numpy arrays
	groups = np.array(groups)
	models = np.array(models)
	betas = np.array(betas)
	sigmas = np.array(sigmas)
	
	# Return
	return(groups, models, betas, sigmas)
	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the FittedISMIP icesheets fitting stage.",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	FittedISMIP_fit_icesheets(args.pipeline_id)
	
	sys.exit()