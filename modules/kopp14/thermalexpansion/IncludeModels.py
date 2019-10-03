import numpy as np
import os

''' IncludeModels.py

This script parses through a directory of models to determine if the output from each
model should be included in the analysis.

Parameters: 
model_dir = Directory of model output. Each model is a subdirectory within this one.
varnames =  Names of the variables of interest. Will be subdirectories of individual 
			model directories.
years = Years of interest.

Return: 
model_list = Vector of model names that are to be included (nmodels)
ZOSTOGA = Global average thermosteric sea-level change (years, nmodels)

Note: The originial code (Matlab version from Kopp14) includes both ZOSTOGA and ZOSGA, 
though it's not apparent to me why these are being treated as the same quantity. ZOSTOGA
is the total thermosteric sea-level change while ZOSGA is total sea-level change including
mass flux from land-ice.

'''

def IncludeModels(model_dir, varnames, years):
	
	# Initialize the model list and data matrix
	model_list = []
	init_zostoga = True
	
	# Loop through available models in model_dir
	for model in os.listdir(model_dir):
		
		# Skip if the folder/file found in this directory is hidden
		if(model.startswith(".")):
			continue
		
		# Initialize "foundvar"
		foundvar = False
		
		# Loop through the variable names
		for this_var in varnames:
			
			# Move onto the next iteration if the variable doesn't exist
			# or if we've already found a viable variable candidate
			vardir = os.path.join(model_dir, model, this_var.lower())
			if(os.path.isdir(vardir) and (not foundvar)):
				
				# Note that we found a viable variable
				foundvar = True
				
				# Note to incorporate this model/variable
				incorporate = True
				
				# Read in the data file
				infile = os.listdir(vardir)[0]
				with open(os.path.join(vardir, infile), 'r') as f:
					try:
						dat = np.loadtxt(f, skiprows=1, delimiter=',')
					except:
						print("Failed to read in data for {0}/{1}".format(vardir, infile))
						incorporate = False
				
				# If the change is too small, ignore this
				good_inds = np.nonzero(~np.isnan(dat[:,1]))[0]
				totalchange = np.sum(np.diff(dat[good_inds[0]:good_inds[-1], 1]))/(dat[good_inds[-1],0] - dat[good_inds[0],0])
				if(totalchange < 2e-4):
					incorporate=False
				
				# If there's no data near year 2000, ignore this
				if(np.min(np.abs(dat[good_inds,0]-2000)) > 2):
					incorporate=False
				
				# If it's a keeper, add it to the list and output data matrix
				if(incorporate):
					model_list.append(model)
					data_to_append = np.interp(years, dat[:,0], dat[:,1], left=np.nan, right=np.nan)
					if(init_zostoga):
						ZOSTOGA = data_to_append
						init_zostoga = False
					else:
						ZOSTOGA = np.vstack((ZOSTOGA, data_to_append))
	
	return(model_list, np.transpose(ZOSTOGA))