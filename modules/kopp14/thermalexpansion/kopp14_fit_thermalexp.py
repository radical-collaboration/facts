import numpy as np
import os
import pickle
import argparse

''' kopp14_fit_thermalexp.py

This script runs the thermal expansion fitting task for the Kopp 2014 workflow. 
There's not much to fit or calibrate.  This task basically aggregates the thermal
expansion quantities into a mean and standard deviation (along with drift correction
if necessary).

Parameters: 
pipeline_id = Unique identifier for the pipeline running this code

Output: Pickle file containing the multi-model mean and standard deviation of global
		thermal expansion. The file also contains the years across which these data are
		available as well as the number of models that contribute to these quantities.

'''

def kopp14_fit_thermalexp(pipeline_id):
	
	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	#rcp_scenario = my_config["rcp_scenario"]
	datayears = my_config["datayears"]
	#targyears = my_config["targyears"]
	#mergeZOSZOSTOGA = my_config["mergeZOSZOSTOGA"]
	#smoothwin = my_config["smoothwin"]
	driftcorr = my_config["driftcorr"]
	#baseyear = my_config["baseyear"]
	#GCMprobscale = my_config["GCMprobscale"]
	
	# Load the ZOSTOGA file
	zostogafile = "{}_ZOSTOGA.pkl".format(pipeline_id)
	try:
		f = open(zostogafile, 'rb')
	except:
		print("Cannot open ZOSTOGA file\n")
	
	# Extract the configuration variables
	my_zostoga = pickle.load(f)
	f.close()

	ZOSTOGA = my_zostoga["ZOSTOGA"]
	histGICrate = my_zostoga["histGICrate"]	
	#CWdrift = my_zostoga["CWdrift"]
	selectyears = my_zostoga["selectyears"]
	
	# Get the mean, std, and counts of models for the thermal expansion component
	# Note: ZOSTOGA mean and std are converted into mm
	ThermExpYears = datayears
	ThermExpMean = np.nanmean(ZOSTOGA, axis=1)*1000
	ThermExpStd = np.nanstd(ZOSTOGA, axis=1)*1000
	def CountNonNAN(x):
		return(len(np.flatnonzero(~np.isnan(x))))
	ThermExpN = np.apply_along_axis(CountNonNAN, axis=1, arr=ZOSTOGA)
	if(driftcorr):
		# NOTE: THIS MAY BE A BUG IN THE OPRIGINAL CODE!!!
		# ThermExpStd has units of mm, but histGICrate has units of meters 
		ThermExpStd = np.sqrt(ThermExpStd**2 + (np.nanstd(histGICrate)*(ThermExpYears-selectyears[0]))**2)  # Consistent with original code
		#ThermExpStd = np.sqrt(ThermExpStd**2 + (np.nanstd(histGICrate*1000)*(ThermExpYears-selectyears[0]))**2)  # Fix the unit mis-match
	
	# Stretch out the thermal expansion metrics for one additional year
	ThermExpMean = np.append(ThermExpMean, (ThermExpMean[-1]+(ThermExpMean[-1]-ThermExpMean[-2])))
	ThermExpStd = np.append(ThermExpStd, (ThermExpStd[-1]+(ThermExpStd[-1]-ThermExpStd[-2])))
	ThermExpYears = np.append(ThermExpYears, (ThermExpYears[-1]+(ThermExpYears[-1]-ThermExpYears[-2])))
	ThermExpN = np.append(ThermExpN, (ThermExpN[-1]))
	
	# Store the thermal expansion variables in a pickle
	output = {'ThermExpMean': ThermExpMean, 'ThermExpStd': ThermExpStd,\
		'ThermExpYears': ThermExpYears,	'ThermExpN': ThermExpN}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the thermal expansion fitting stage for the Kopp14 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the fitting process on the file specified from the command line argument
	kopp14_fit_thermalexp(args.pipeline_id)
	
	exit()