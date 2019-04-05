import numpy as np
import pickle as p
import os

''' kopp14_preprocess_icesheets.py

This script runs the ice sheet pre-processing task for the Kopp 2014 workflow. Currently,
the data are hard-coded in this script, from which a pickled file is generated. This
could easily be altered to take filenames as parameters from which these data are read. 
This would be useful for testing multiple sets of these data (likely with different
assumptions) in batch.

Parameters: NONE

Output: Pickled file containing barates, lastdecadegt, and aris2090 variables within a
dictionary vairable.

'''

if __name__ == '__main__':

	# Fill the rates data matricies
	barates = np.array([[0.8, 1.0, 1.2, 2.4, 5.8],\
						[0.2, 3.0, 0.3, 1.5, 11.8],\
						[-1.9, 2.8, -1.5, 0.2, 10.2]])
	lastdecadegt = np.array([-211, -85-29, 26])
	aris2090 = np.array([[[0.07, 0.12, 0.21], [-0.06, 0.04, 0.12]],\
						[[0.04, 0.08, 0.13], [-0.04, 0.05, 0.13]],\
						[[0.04, 0.08, 0.13], [-0.04, 0.05, 0.13]],\
						[[0.04, 0.06, 0.10], [-0.03, 0.05, 0.14]]]) * 1000
					
	# Fill the correlation data matricies
	bacorris = np.array([[1.0, 0.7, -0.2], [0.7, 1.0, -0.2], [-0.2, -0.2, 1]])
	arcorris = np.array([[1.0, 0.5], [0.5, 1.0]])
	#bacorris = np.eye(3)
	#arcorris = np.eye(2)

	# Collate rates data into a single dictionary
	data_is = {'barates': barates, 'lastdecadegt': lastdecadegt, 'aris2090': aris2090[0,:,:]}

	# Collate the correlation data into a single dictionary
	corr_is = {'bacorris': bacorris, 'arcorris': arcorris}
	
	# Define the data directory
	outdir = os.path.join(os.path.dirname(__file__), 'data')

	# Write the rates data to a pickle file
	ratefile = open(os.path.join(outdir, "kopp14_icesheets_rates.pkl"), 'wb')
	p.dump(data_is, ratefile)
	ratefile.close()

	# Write the correlation data to a pickle file
	corrfile = open(os.path.join(outdir, "kopp14_icesheets_corr.pkl"), 'wb')
	p.dump(corr_is, corrfile)
	corrfile.close() 

	# Done!
	exit()