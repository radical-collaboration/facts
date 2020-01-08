import numpy as np
import os
import pickle
import argparse
import re

''' kopp14_fit_oceandynamics.py

This runs the fitting stage for the ocean dynamics component of the Kopp14 workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

'''

def kopp14_fit_oceandynamics(pipeline_id):

	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open configuration file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	datayears = my_config["datayears"]
	mergeZOSZOSTOGA = my_config["mergeZOSZOSTOGA"]
	smoothwin = my_config["smoothwin"]
	driftcorr = my_config["driftcorr"]
	baseyear = my_config["baseyear"]
	
	# Load the ZOSTOGA file
	zostogafile = "{}_ZOSTOGA.pkl".format(pipeline_id)
	try:
		f = open(zostogafile, 'rb')
	except:
		print("Cannot open ZOSTOGA file\n")
	
	# Extract the ZOSTOGA variables
	my_zostoga = pickle.load(f)
	f.close()

	sZOSTOGA = my_zostoga["sZOSTOGA"]
	histGICrate = my_zostoga["histGICrate"]	
	zostoga_modellist = my_zostoga['zostoga_modellist']
	selectyears = my_zostoga["selectyears"]
	
	# Load the ZOS file
	zosfile = "{}_ZOS.pkl".format(pipeline_id)
	try:
		f = open(zosfile, 'rb')
	except:
		print("Cannot open ZOS file\n")
	
	# Extract the ZOS variables
	my_zos = pickle.load(f)
	f.close()

	sZOS = my_zos["sZOS"]
	sZOSTOGAadj = my_zos["sZOSTOGAadj"]
	focus_site_lats = my_zos["focus_site_lats"]
	focus_site_ids = my_zos["focus_site_ids"]
	comb_modellist = my_zos["comb_modellist"]
	
	
	#------------- Begin Thermal Expansion -----------------------------------------------
	
	# Get the mean, std, and counts of models for the thermal expansion component
	# Note: ZOSTOGA mean and std are converted into mm
	ThermExpYears = datayears
	ThermExpMean = np.nanmean(sZOSTOGA, axis=1)*1000
	ThermExpStd = np.nanstd(sZOSTOGA, axis=1)*1000
	def CountNonNAN(x):
		return(len(np.flatnonzero(~np.isnan(x))))
	ThermExpN = np.apply_along_axis(CountNonNAN, axis=1, arr=sZOSTOGA)
	if(driftcorr):
		# NOTE: THIS MAY BE A BUG IN THE ORIGINAL CODE!!!
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
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_thermalexp_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	
	#-------------------- Begin Ocean Dynamics -------------------------------------------
	
	# Pick out indices for year between 2000 and 2100 (consistent with K14) and trim data
	year_idx = np.flatnonzero(np.logical_and(datayears > 2000, datayears < 2100))
	sZOS = sZOS[year_idx,:,:]
	
	# Calculate and remove "extremeness" as models whose year 2099 value is over 
	# 10x the median across models in year 2099
	ext_num = np.abs(sZOS[-1,:,:])
	ext_denom = np.nanmedian(np.abs(sZOS[-1,:,:]), axis=0)
	extremeness = ext_num/ext_denom
	with np.errstate(invalid='ignore'):
		model_idx = extremeness < 10   # Wrapped in np.errstate call to surpress warning of 'nan' in less-than test
	nan_mask = np.where(model_idx, 1, np.nan)
	sZOS = sZOS * nan_mask[np.newaxis,:,:]	
	
	# Calculate and remove "extremeness" as models that in year 2099 have values
	# (less the mean across models) are greater than 3 standard deviations across models
	ext_num = sZOS[-1,:,:] - np.nanmean(sZOS[-1,:,:], axis=0)
	ext_denom = np.nanstd(sZOS[-1,:,:], axis=0)
	extremeness = np.abs(ext_num/ext_denom)
	std_limit = np.where(ext_denom > 0.2, 1, 0)
	with np.errstate(invalid='ignore'):
		model_idx = extremeness > 3   # Wrapped in np.errstate call to surpress warning of 'nan' in greater-than test
	nan_mask = np.where(model_idx * std_limit == 1, np.nan, 1)
	sZOS = sZOS * nan_mask[np.newaxis,:,:]
	
	# If the latitude for the sites of interest is greater than 50, remove the "miroc"
	# and "giss" models
	lat_limit = np.where(focus_site_lats > 50.0, 1, 0)
	model_idx = np.array([bool(re.search(r"miroc|giss", x)) for x in comb_modellist])
	nan_mask = np.where(model_idx[:,np.newaxis] * lat_limit == 1, np.nan, 1)
	sZOS = sZOS * nan_mask[np.newaxis,:,:]

	# Calculate the OD mean, std, and N
	OceanDynMean = np.nanmean(sZOS, axis=1) * 1000
	OceanDynStd = np.nanstd(sZOS, axis=1) * 1000
	OceanDynN = np.nansum(~np.isnan(sZOS), axis=1)
	
	# Trim sZOSTOGAadj to same year range as sZOS
	sZOSTOGAadj = sZOSTOGAadj[year_idx,:]
	
	# Calculate the correlation of ZOS with thermal expansion
	zos_demean = sZOS - np.nanmean(sZOS, axis=1)[:,np.newaxis,:]
	zostoga_demean = sZOSTOGAadj - np.nanmean(sZOSTOGAadj, axis=1)[:,np.newaxis]
	corr_num = np.nansum(zos_demean * zostoga_demean[:,:,np.newaxis], axis=1)
	corr_denom = np.sqrt(np.nansum(zos_demean**2, axis=1) * np.nansum(zostoga_demean**2, axis=1)[:,np.newaxis])
	OceanDynTECorr = corr_num / corr_denom	
	
	# Extend the OceanDyn* variables to 2100
	OceanDynMean = np.append(OceanDynMean, OceanDynMean[-1,:] + (OceanDynMean[-1,:] - OceanDynMean[-2,:])[np.newaxis,:], axis=0)
	OceanDynStd = np.append(OceanDynStd, OceanDynStd[-1,:] + (OceanDynStd[-1,:] - OceanDynStd[-2,:])[np.newaxis,:], axis=0)
	OceanDynN = np.append(OceanDynN, (OceanDynN[-1,:])[np.newaxis,:], axis=0)
	OceanDynTECorr = np.append(OceanDynTECorr, OceanDynTECorr[-1,:] + (OceanDynTECorr[-1,:] - OceanDynTECorr[-2,:])[np.newaxis,:], axis=0)
	
	# Ensure correlation remains within [-1,1]
	OceanDynTECorr = np.maximum(-1.0, np.minimum(1.0, OceanDynTECorr))
	
	# Define the years over which Ocean Dynamics mean, std, N, and TECorr are defined
	OceanDynYears = np.append(datayears[year_idx], datayears[year_idx[-1]] + 1)
	
	# Store the ocean dynamics variables in a pickle
	output = {'OceanDynMean': OceanDynMean, 'OceanDynStd': OceanDynStd,\
		'OceanDynYears': OceanDynYears,	'OceanDynN': OceanDynN, 'OceanDynTECorr': OceanDynTECorr}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_oceandynamics_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	
if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the Kopp14 ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Kopp14 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the user defined RCP scenario
	kopp14_fit_oceandynamics(args.pipeline_id)
	
	# Done
	exit()