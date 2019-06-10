import numpy as np
import os
import pickle

''' kopp14_fit_thermalexp.py

This script runs the thermal expansion fitting task for the Kopp 2014 workflow. 
There's not much to fit or calibrate.  This task basically aggregates the thermal
expansion quantities into a mean and standard deviation (along with drift correction
if necessary).

Parameters: 
configfile = Configuration pickle produced by kopp14_preprocess_thermalexp.py
zostogafile = ZOSTOGA pickle produced b kopp14_preprocess_thermalexp.py

Output: Pickle file containing the multi-model mean and standard deviation of global
		thermal expansion. The file also contains the years across which these data are
		available as well as the number of models that contribute to these quantities.

'''

def kopp14_fit_thermalexp(configfile, zostogafile):
	
	# Load the configuration file
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	#scenario = my_config["scenario"]
	datayears = my_config["datayears"]
	#targyears = my_config["targyears"]
	#mergeZOSZOSTOGA = my_config["mergeZOSZOSTOGA"]
	#smoothwin = my_config["smoothwin"]
	driftcorr = my_config["driftcorr"]
	#baseyear = my_config["baseyear"]
	#GCMprobscale = my_config["GCMprobscale"]
	
	# Load the ZOSTOGA file
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
	#selectyears = my_zostoga["selectyears"]
	
	# Get the mean, std, and counts of models for the thermal expansion component
	ThermExpYears = datayears
	ThermExpMean = np.nanmean(ZOSTOGA, axis=1)*1000
	ThermExpStd = np.nanstd(ZOSTOGA, axis=1)*1000
	def CountNonNAN(x):
		return(len(np.flatnonzero(~np.isnan(x))))
	ThermExpN = np.apply_along_axis(CountNonNAN, axis=1, arr=ZOSTOGA)
	if(driftcorr):
		ThermExpStd = np.sqrt(ThermExpStd**2 + (np.nanstd(histGICrate)*(ThermExpYears-selectyears[0]))**2)
	
	# Stretch out the thermal expansion metrics for one additional year
	ThermExpMean = np.append(ThermExpMean, (ThermExpMean[-1]+(ThermExpMean[-1]-ThermExpMean[-2])))
	ThermExpStd = np.append(ThermExpStd, (ThermExpStd[-1]+(ThermExpStd[-1]-ThermExpStd[-2])))
	ThermExpYears = np.append(ThermExpYears, (ThermExpYears[-1]+(ThermExpYears[-1]-ThermExpYears[-2])))
	ThermExpN = np.append(ThermExpN, (ThermExpN[-1]))
	
	# Store the thermal expansion variables in a pickle
	output = {'ThermExpMean': ThermExpMean, 'ThermExpStd': ThermExpStd,\
		'ThermExpYears': ThermExpYears,	'ThermExpN': ThermExpN}
	outfile = open(os.path.join(os.path.dirname(__file__), "data", "kopp14_thermalexp_fit.pkl"), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

if __name__ == '__main__':
	
	# Run the fitting process on the file specified from the command line argument
	try:
		config_file = sys.argv[1]
	except:
		config_file = os.path.join(os.path.dirname(__file__), "data", "kopp14_thermalexp_config.pkl")
		
	try:
		zostogafile = sys.argv[2]
	except:
		zostogafile = os.path.join(os.path.dirname(__file__), "data", "kopp14_thermalexp_ZOSTOGA.pkl")
	finally:
		kopp14_fit_thermalexp(config_file, zostogafile)
	
	exit()