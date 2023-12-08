import numpy as np
import os
import sys
import pickle
import argparse
import re

''' tlm_fit_oceandynamics.py

This runs the fitting stage for the ocean dynamics and thermal expansion
component of the TLM workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

'''

def tlm_fit_oceandynamics(pipeline_id):

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
	no_correlation = my_config["no_correlation"]
	maxDOF = my_config["maxDOF"]

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
	zostoga_modellist = my_zostoga['zostoga_modellist']

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
	zos_modellist = my_zos["zos_modellist"]

	# Subset of data years
	year_idx = np.flatnonzero(np.logical_and(datayears > 2000, datayears <= 2300))
	year_2099_idx = np.flatnonzero(datayears[year_idx] == 2099)

	# Function to replace locations/times where number of models is below a
	# threshold with an interpolation/extrapolation.
	def MyInterp(x, xp, fp, extrap=False, right=np.nan, left=np.nan, no_neg_rate=False):
		sort_idx = np.argsort(xp)
		temp_xp = xp[sort_idx]
		temp_fp = fp[sort_idx]
		delta_idx = np.flatnonzero(np.logical_and(temp_xp > 2060, temp_xp <= 2100))
		extrap_delta = np.mean(np.diff(temp_fp[delta_idx])/np.diff(temp_xp[delta_idx]))
		if extrap_delta < 0.0 and no_neg_rate:
			extrap_delta = 0.0
		def pointwise(x):
			if x < temp_xp[0]:
				return(temp_fp[0] + (x-temp_xp[0])*extrap_delta)
			if x > temp_xp[-1]:
				return(temp_fp[-1] + (x-temp_xp[-1])*extrap_delta)
			else:
				return(np.interp(x, temp_xp, temp_fp))
		if extrap:
			return(np.array(list(map(pointwise, x))))
		else:
			return(np.interp(x, xp, fp, left=left, right=right))

	# Minimum number of models needed for mean/sd/corr calculations
	# To Do: This should eventually be a commandline option
	min_n_models = 5

	# Apply extrapolation?
	# To Do: This should eventually be a commandline option
	do_extrap = True

	#------------- Begin Thermal Expansion -----------------------------------------------

	# Get the mean, std, and counts of models for the thermal expansion component
	# Note: ZOSTOGA mean and std are converted into mm
	ThermExpYears = datayears[year_idx]
	ThermExpMean = np.nanmean(sZOSTOGA[year_idx,:], axis=1)*1000
	ThermExpStd = np.nanstd(sZOSTOGA[year_idx,:], axis=1)*1000
	def CountNonNAN(x):
		return(len(np.flatnonzero(~np.isnan(x))))
	ThermExpN = np.apply_along_axis(CountNonNAN, axis=1, arr=sZOSTOGA[year_idx,:])
	#if(driftcorr):
		# NOTE: THIS MAY BE A BUG IN THE ORIGINAL CODE!!!
		# ThermExpStd has units of mm, but histGICrate has units of meters
		#ThermExpStd = np.sqrt(ThermExpStd**2 + (np.nanstd(histGICrate)*(ThermExpYears-selectyears[0]))**2)  # Consistent with original code
		#ThermExpStd = np.sqrt(ThermExpStd**2 + (np.nanstd(histGICrate*1000)*(ThermExpYears-selectyears[0]))**2)  # Fix the unit mis-match

	# For consistency with original matlab code, apply linear extrapolation to year 2100
	#ThermExpMean[year_2099_idx+1] = ThermExpMean[year_2099_idx]+(ThermExpMean[year_2099_idx]-ThermExpMean[year_2099_idx-1])
	#ThermExpStd[year_2099_idx+1] = ThermExpStd[year_2099_idx]+(ThermExpStd[year_2099_idx]-ThermExpStd[year_2099_idx-1])
	#ThermExpN[year_2099_idx+1] = ThermExpN[-1]

	# Determine the DOF
	ThermExpDOF = np.copy(ThermExpN)

	# Interpolate/extrapolate values as needed
	temp_good_idx = np.flatnonzero(ThermExpN >= min_n_models)
	if len(temp_good_idx) > 0:
		temp_replace_idx = np.flatnonzero(ThermExpN <= min_n_models)
		ThermExpDOF[temp_replace_idx] = ThermExpN[temp_good_idx[-1]]
		ThermExpStd[temp_replace_idx] = MyInterp(ThermExpYears[temp_replace_idx], ThermExpYears[temp_good_idx], ThermExpStd[temp_good_idx], extrap=do_extrap, right=ThermExpStd[temp_good_idx[-1]], no_neg_rate=True)
		ThermExpMean[temp_replace_idx] = MyInterp(ThermExpYears[temp_replace_idx], ThermExpYears[temp_good_idx], ThermExpMean[temp_good_idx], extrap=do_extrap, right=ThermExpMean[temp_good_idx[-1]])

	else:
		ThermExpStd[:] = 0.0
		ThermExpMean[:] = 0.0

	# Ensure valid DOF values
	ThermExpDOF -= 1
	ThermExpDOF[ThermExpDOF > maxDOF] = maxDOF
	ThermExpDOF[ThermExpDOF < 1] = 1

	# Store the thermal expansion variables in a pickle
	output = {'ThermExpMean': ThermExpMean, 'ThermExpStd': ThermExpStd,\
		'ThermExpYears': ThermExpYears,	'ThermExpN': ThermExpN, 'ThermExpDOF': ThermExpDOF}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_thermalexp_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=4)
	outfile.close()


	#-------------------- Begin Ocean Dynamics -------------------------------------------

	# Pick out indices for year between 2000 and 2300 (consistent with K14) and trim data
	sZOS = sZOS[year_idx,:,:]

	# Determine which locations have enough models initially to warrant further
	# checks on extremeness. Points with less than 7 models available are
	# permitted to bypass the extremeness checks.
	extremeness_model_check = np.sum(~np.isnan(sZOS[year_2099_idx,:,:]), axis=1) < 7

	# For points that have enough pre-extremeness check models, calculate and
	# remove "extremeness" as models whose year 2099 value is over 10x the
	# median across models in year 2099.
	ext_num = np.abs(sZOS[year_2099_idx,:,:])
	ext_denom = np.nanmedian(np.abs(sZOS[year_2099_idx,:,:]), axis=1)
	extremeness = ext_num/ext_denom
	with np.errstate(invalid='ignore'):
		model_idx = extremeness < 10   # Wrapped in np.errstate call to surpress warning of 'nan' in less-than test
	nan_mask = np.where(np.logical_or(model_idx, extremeness_model_check), 1.0, np.nan)
	sZOS = sZOS * nan_mask

	# For points that have enough pre-extremeness check models, calculate and
	# remove "extremeness" as models that in year 2099 have values (less the
	# mean across models) are greater than 3 standard deviations across models
	ext_num = sZOS[year_2099_idx,:,:] - np.nanmean(sZOS[year_2099_idx,:,:], axis=1)
	ext_denom = np.nanstd(sZOS[year_2099_idx,:,:], axis=1)
	extremeness = np.abs(ext_num/ext_denom)
	std_limit = np.where(ext_denom < 0.2, 1, 0)
	with np.errstate(invalid='ignore'):
		model_idx = extremeness < 3   # Wrapped in np.errstate call to surpress warning of 'nan' in greater-than test
	nan_mask = np.where(np.logical_or(model_idx, np.logical_or(extremeness_model_check,std_limit)), 1.0, np.nan)
	sZOS = sZOS * nan_mask

	# Calculate the OD mean, std, and N
	OceanDynMean = np.nanmean(sZOS, axis=1) * 1000.0
	OceanDynStd = np.nanstd(sZOS, axis=1) * 1000.0
	OceanDynN = np.nansum(~np.isnan(sZOS), axis=1)
	OceanDynDOF = np.copy(OceanDynN)

	# Any points that have 0 models available for the calculation should be have
	# the mean set to NAN
	OceanDynMean[OceanDynN == 0] = np.nan

	# Trim sZOSTOGAadj to same year range as sZOS
	sZOSTOGAadj = sZOSTOGAadj[year_idx,:]

	# Calculate the correlation of ZOS with thermal expansion if needed
	# Note: Correlation returns 'np.nan' if 'corr_denom' == 0
	if no_correlation:
		OceanDynTECorr = 0.0
	else:
		zos_demean = sZOS - np.nanmean(sZOS, axis=1)[:,np.newaxis,:]
		zostoga_demean = sZOSTOGAadj - np.nanmean(sZOSTOGAadj, axis=1)[:,np.newaxis]
		corr_num = np.nansum(zos_demean * zostoga_demean[:,:,np.newaxis], axis=1)
		corr_denom = np.sqrt(np.nansum(zos_demean**2, axis=1) * np.nansum(zostoga_demean**2, axis=1)[:,np.newaxis])
		OceanDynTECorr = corr_num / corr_denom

	# For consistency with original matlab code, apply the year 2100 extension
	#OceanDynMean[year_2099_idx+1,:] = OceanDynMean[year_2099_idx,:] + (OceanDynMean[year_2099_idx,:] - OceanDynMean[year_2099_idx-1,:])
	#OceanDynStd[year_2099_idx+1,:] = OceanDynStd[year_2099_idx,:] + (OceanDynStd[year_2099_idx,:] - OceanDynStd[year_2099_idx-1,:])
	#OceanDynN[year_2099_idx+1,:] = OceanDynN[year_2099_idx+1,:]
	#if not no_correlation:
		#OceanDynTECorr[year_2099_idx+1,:] = OceanDynTECorr[year_2099_idx,:] + (OceanDynTECorr[year_2099_idx,:] - OceanDynTECorr[year_2099_idx-1,:])

	# Define the years over which Ocean Dynamics mean, std, N, and TECorr are defined
	OceanDynYears = datayears[year_idx]

	# Set any instance where correlation with TE is NAN to zero
	OceanDynTECorr[np.isnan(OceanDynTECorr)] = 0.0

	# Interpolate/extrapolate values as needed
	for i in np.arange(len(focus_site_ids)):
		temp_good_idx = np.flatnonzero(OceanDynN[:,i] >= min_n_models)
		if len(temp_good_idx) > 0:
			temp_replace_idx = np.flatnonzero(OceanDynN[:,i] <= min_n_models)
			OceanDynDOF[temp_replace_idx,i] = OceanDynN[temp_good_idx[-1],i]
			OceanDynTECorr[temp_replace_idx,i] = MyInterp(OceanDynYears[temp_replace_idx], OceanDynYears[temp_good_idx], OceanDynTECorr[temp_good_idx,i], right=OceanDynTECorr[temp_good_idx[-1],i])
			OceanDynStd[temp_replace_idx,i] = MyInterp(OceanDynYears[temp_replace_idx], OceanDynYears[temp_good_idx], OceanDynStd[temp_good_idx,i], extrap=do_extrap, right=OceanDynStd[temp_good_idx[-1],i], no_neg_rate=True)
			OceanDynMean[temp_replace_idx,i] = MyInterp(OceanDynYears[temp_replace_idx], OceanDynYears[temp_good_idx], OceanDynMean[temp_good_idx,i], extrap=do_extrap, right=OceanDynMean[temp_good_idx[-1],i])

		else:
			OceanDynStd[:,i] = 0.0
			OceanDynMean[:,i] = np.nan	#0.0
			OceanDynTECorr[:,i] = 0.0

	# Ensure correlation remains within [-1,1]
	OceanDynTECorr = np.maximum(-1.0, np.minimum(1.0, OceanDynTECorr))

	# Ensure valid DOF values
	OceanDynDOF -= 1
	OceanDynDOF[OceanDynDOF > maxDOF] = maxDOF
	OceanDynDOF[OceanDynDOF < 1] = 1

	# Store the ocean dynamics variables in a pickle
	output = {'OceanDynMean': OceanDynMean, 'OceanDynStd': OceanDynStd, 'OceanDynDOF': OceanDynDOF, \
		'OceanDynYears': OceanDynYears,	'OceanDynN': OceanDynN, 'OceanDynTECorr': OceanDynTECorr}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_oceandynamics_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile,protocol=4)
	outfile.close()

if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the TLM ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Define the names of the intermediate files this run will generate
	outdir = os.path.dirname(__file__)
	odfitfile = os.path.join(outdir, "{}_oceandynamics_fit.pkl".format(args.pipeline_id))
	tefitfile = os.path.join(outdir, "{}_thermalexp_fit.pkl".format(args.pipeline_id))

	# Run the OD fitting stage if intermediate files are not present
	if os.path.isfile(tefitfile) and os.path.isfile(odfitfile):
		print("{}, {} found, skipping OD fitting".format(tefitfile, odfitfile))
	else:
		tlm_fit_oceandynamics(args.pipeline_id)


	# Done
	sys.exit()
