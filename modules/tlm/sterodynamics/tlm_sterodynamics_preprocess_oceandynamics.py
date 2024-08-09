import numpy as np
import pickle
import os
import sys
import re
import argparse
import fnmatch
from IncludeCMIP6Models import IncludeCMIP6Models
from IncludeCMIP6ZOSModels import *
from SmoothZOSTOGA import SmoothZOSTOGA
#from DriftCorr import DriftCorr
from read_locationfile import ReadLocationFile
from Smooth import Smooth

''' tlm_preprocess_oceandynamics.py

This runs the preprocessing stage for the ocean dynamics component of the IPCC AR6
workflow.

Parameters:
ssp_scenario = SSP scenario of interest
modeldir = Directory that contains the ZOS and ZOSTOGA CMIP6 model data
driftcorr = Apply the drift correction?
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code


'''

def FindInputModels(tasdir, zosdir, scenario):

	# Acceptable SSP scenarios
	ssp_scenarios = ['ssp585', 'ssp370', 'ssp245', 'ssp126', 'ssp119']

	# Initialize models and scenarios to include
	include_models = []
	include_scenarios = []

	# Test the provided scenario
	scenario_test = re.search("^tlim(\d*\.?\d+)win(\d*\.?\d+)$", scenario)
	if(scenario_test):

		# This is a temperature limit, so extract the limit from the scenario string
		temp_target = float(scenario_test.group(1))
		temp_target_window = float(scenario_test.group(2))

		# Produce a list of models and scenarios that match the criteria
		(include_models, include_scenarios) = tas_limit_filter(tasdir, temp_target, temp_target_window)

	# This scenario has no need for identifying models with TAS variable overlap.
	# Simply return a list of zoz models available for the selected scenario
	elif(scenario in ssp_scenarios):

		# Loop through the list of models available in the tas directory
		for this_model in os.listdir(zosdir):

			# Skip hidden directories and files
			if(re.search("^\.", this_model)):
				continue

			# Find the appropriate zos file for this model and scenario
			model_has_scenario = fnmatch.filter(os.listdir(os.path.join(zosdir, this_model)), "*{}*".format(scenario))
			if(not model_has_scenario):
				continue
			zos_filename = model_has_scenario[0]

			# If the file exists, append this model to the model list
			if(os.path.isfile(os.path.join(zosdir, this_model, zos_filename))):
				include_models.append(this_model)
				include_scenarios.append(scenario)

	else:

		# This is an invalid scenario
		raise Exception("Invalid scenario definition: {}".format(scenario))

	return(include_models, include_scenarios)



def tas_limit_filter(tasdir, temp_target, temp_target_window=0.25, ref_syear=1850, ref_eyear=1900, extrap_eyear=2099):

	# Initialize a running list of models and scenarios to include in the analysis
	include_models = []
	include_scenarios = []

	# Get a list of models available from the subdirectories available in the parent model directory
	modeldirs = os.listdir(tasdir)

	# Loop through the models
	for this_modeldir in modeldirs:

		# Skip any hidden directories or files
		if(re.search("^\.", this_modeldir)):
			continue

		# Open the historical file if available
		hist_filename = fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*historical*")[0]
		try:
			nc_hist = Dataset(os.path.join(tasdir, this_modeldir, hist_filename), 'r')
		except:
			continue

		# Get the years out of the historical file
		#hist_years = nc_hist.variables['year'][:]
		nctime = nc_hist.variables["time"]
		temp_nctime = cftime.num2date(nctime, nctime.units, nctime.calendar)
		hist_years = np.array([int(x.strftime("%Y")) for x in temp_nctime])

		# Determine which indices are in our reference average window
		ref_idx = np.flatnonzero(np.logical_and(hist_years >= ref_syear, hist_years <= ref_eyear))

		# Extract the temperature data from the reference period and produce the average
		ref_tas = np.mean(nc_hist.variables['tas'][ref_idx])

		# Close the historical netCDF file
		nc_hist.close()

		# Loop through all the remaining files in this model directory
		for this_filename in fnmatch.filter(os.listdir(os.path.join(tasdir, this_modeldir)), "*ssp*"):

			# Open the netCDF file if available
			this_file = os.path.join(tasdir, this_modeldir, this_filename)
			nc = Dataset(this_file, 'r')

			# Get the years available in this projection
			#proj_years = nc.variables['year'][:]
			nctime = nc.variables["time"]
			temp_nctime = cftime.num2date(nctime, nctime.units, nctime.calendar)
			proj_years = np.array([int(x.strftime("%Y")) for x in temp_nctime])

			# Determine the indices needed to calculate the rate of the 19-yr average
			stop_extrap_idx = np.flatnonzero(np.logical_and(proj_years >= (extrap_eyear-19), proj_years <= extrap_eyear))
			start_extrap_idx = np.flatnonzero(np.logical_and(proj_years >= (extrap_eyear-20), proj_years <= (extrap_eyear-1)))

			# Calculate the means and take the difference to get the rate
			start_extrap = np.mean(nc.variables['tas'][start_extrap_idx])
			stop_extrap = np.mean(nc.variables['tas'][stop_extrap_idx])
			tas_rate = stop_extrap - start_extrap

			# Extrapolate that to get a 19-yr average centered on 2100 and subtract the reference
			tas_mean_2100 = (stop_extrap + (2100 - extrap_eyear + 9) * tas_rate) - ref_tas

			'''
			# Get the annual mean tas values
			tas_vals = nc.variables['tas']
			tas_2100 = (tas_vals[-1] - tas_vals[-2]) + tas_vals[-1]
			tas_vals = np.append(tas_vals, tas_2100)
			proj_years = np.append(proj_years, 2100)
			'''


			# If tas_mean_2100 is within 0.25 deg of the temperature target, set this model as a keeper
			if(tas_mean_2100 >= temp_target - temp_target_window and tas_mean_2100 <= temp_target + temp_target_window):

				# Add the model to the list
				include_models.append(this_modeldir)

				# Add the scenario to the list
				this_scenario = re.search(r"(ssp\d{3})", this_filename).group(1)
				include_scenarios.append(this_scenario)

	return(include_models, include_scenarios)



def tlm_preprocess_oceandynamics(scenario, modeldir, driftcorr, no_correlation, pyear_start, pyear_end, pyear_step, locationfilename, baseyear, pipeline_id):

	# Define variables
	datayears = np.arange(1861,2301)
	targyears = np.arange(pyear_start,pyear_end+1,pyear_step)
	smoothwin = 19
	GCMprobscale = 0.833
	maxDOF = np.iinfo(np.int32).max

	# Define the list of input models and scenarios
	tasdir = os.path.join(modeldir, "tas")
	zos_modeldir = os.path.join(modeldir, "zos")
	zostoga_modeldir = os.path.join(modeldir, "zostoga")
	(include_models, include_scenarios) = FindInputModels(tasdir, zos_modeldir, scenario)

	if(not include_models):
		raise Exception("No models found for this scenario or temperature target - {}".format(scenario))

	# Turn off correlation if this is a temperature target run
	#if(re.search("^tlim", scenario)):
	#	no_correlation = True

	# Turn off merging ZOS and ZOSTOGA if no_correlation
	#if no_correlation:
	#	mergeZOSZOSTOGA = 0
	#else:
	#	mergeZOSZOSTOGA = 1

	# Merging is done in the postprocessing stage automatically.
	mergeZOSZOSTOGA = 0

	# Read in the ZOSTOGA data
	(zostoga_modellist, zostoga_scenariolist, ZOSTOGA) = IncludeCMIP6Models(zostoga_modeldir,'zostoga', datayears, include_models, include_scenarios)

	# Center, suture, and smooth ZOSTOGA
	sZOSTOGA = np.nan * ZOSTOGA
	for i in np.arange(0,ZOSTOGA.shape[1]):
		(ZOSTOGA[:,i], sZOSTOGA[:,i]) = SmoothZOSTOGA(ZOSTOGA[:,i], datayears, baseyear, smoothwin)

	# Store the configuration in a pickle
	output = {'datayears': datayears, 'scenario': scenario, \
		'targyears': targyears,	'mergeZOSZOSTOGA': mergeZOSZOSTOGA,\
		'smoothwin': smoothwin, 'driftcorr': driftcorr, 'baseyear': baseyear,\
		'GCMprobscale': GCMprobscale, 'maxDOF': maxDOF, 'no_correlation': no_correlation}

	# Write the configuration to a file
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_config.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=4)
	outfile.close()

	# Store the ZOSTOGA variables in a pickle
	output = {'sZOSTOGA': sZOSTOGA, 'zostoga_modellist': zostoga_modellist, \
				'zostoga_scenariolist': zostoga_scenariolist}

	# Write the ZOSTOGA variables to a file
	outfile = open(os.path.join(outdir, "{}_ZOSTOGA.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=4)
	outfile.close()


	#------------ Begin Ocean Dynamics ---------------------------------------------------

	# Load the site locations
	locationfile = os.path.join(os.path.dirname(__file__), locationfilename)
	(_, focus_site_ids, focus_site_lats, focus_site_lons) = ReadLocationFile(locationfile)

	# Load the ZOS data
	(zos_modellist, zos_scenariolist, ZOS_raw) = IncludeCMIP6ZOSModels(zos_modeldir, "zos", datayears, include_models, include_scenarios, focus_site_lats, focus_site_lons)

	# Find the overlap between ZOS and ZOSTOGA
	comb_modellist, zostoga_model_idx, zos_model_idx = np.intersect1d(zostoga_modellist, zos_modellist, return_indices=True)

	'''
	#NOTE: POTENTIAL BUG IN ORIGINAL CODE
	#The original code uses the smoothed ZOSTOGA data as the raw ZOSTOGA values. This in
	#turn smooths the ZOSTOGA data over the 'smoothwin' period twice. To replicate the
	#bug, set 'ZOSTOGAadj' to a subset of 'sZOSTOGA' instead of 'ZOSTOGA'.
	'''
	# If no_correlation, do not subset the models to overlap
	if no_correlation:
		ZOSTOGAadj = sZOSTOGA # Replicate potential bug
		#ZOSTOGAadj = ZOSTOGA  # Fix for potential bug
	else:
		ZOSTOGAadj = sZOSTOGA[:,zostoga_model_idx]  # Replicate potential bug
		#ZOSTOGAadj = ZOSTOGA[:,zostoga_model_idx]  # Fix for potential bug
		ZOS_raw = ZOS_raw[:, zos_model_idx, :]

	# Should we merge ZOSTOGA and ZOS?
	if(mergeZOSZOSTOGA):
		ZOS = ZOS_raw + ZOSTOGAadj[:,:,np.newaxis]
	else:
		ZOS = ZOS_raw

	# Smooth ZOS and ZOSTOGA over 19 year smoothing window
	def nanSmooth(x, w):
		idx = np.flatnonzero(~np.isnan(x))
		temp = x
		if len(idx) > 0:
			temp[idx] = Smooth(x[idx], w)
		return(temp)

	sZOS = np.apply_along_axis(nanSmooth, axis=0, arr=ZOS, w=smoothwin)
	sZOSTOGAadj = np.apply_along_axis(nanSmooth, axis=0, arr=ZOSTOGAadj, w=smoothwin)

	# Center the smoothed ZOS/ZOSTOGAadj to the baseyear
	baseyear_idx = np.flatnonzero(datayears == baseyear)
	sZOS = np.apply_along_axis(lambda z, idx: z - z[idx], axis=0, arr=sZOS, idx=baseyear_idx)
	sZOSTOGAadj = np.apply_along_axis(lambda z, idx: z - z[idx], axis=0, arr=sZOSTOGAadj, idx=baseyear_idx)

	# Store the ZOS variable in a pickle
	output = {'sZOS': sZOS, 'zos_modellist': zos_modellist, 'zos_scenariolist': zos_scenariolist, 'datayears': datayears, \
				'focus_site_ids': focus_site_ids, 'focus_site_lats': focus_site_lats, \
				'focus_site_lons': focus_site_lons, 'sZOSTOGAadj': sZOSTOGAadj, 'comb_modellist': comb_modellist}

	# Write the ZOS variables to a file
	outfile = open(os.path.join(outdir, "{}_ZOS.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile, protocol=4)
	outfile.close()



if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the TLM ocean dynamics workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)", default='ssp585')
	parser.add_argument('--scenario_dsl', help="SSP scenario to use for correlation of thermal expansion and dynamic sea level, if not the same as scenario", default='', choices=['','ssp119','ssp126','ssp245','ssp370','ssp585'])
	parser.add_argument('--model_dir', help="Directory containing ZOS/ZOSTOGA CMIP6 GCM output",\
	default=os.path.join(os.path.dirname(__file__), "cmip6"))

	parser.add_argument('--no_drift_corr', help="Do not apply the drift correction", type=int, choices=[0,1], default=0)

	parser.add_argument('--no_correlation', help="Do not apply the correlation between ZOS and ZOSTOGA fields", type=int, choices=[0,1], default=0)

	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)

	parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")

	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2000)

	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()

	# Define the names of the intermediate files this run will generate
	outdir = os.path.dirname(__file__)
	configfile = os.path.join(outdir, "{}_config.pkl".format(args.pipeline_id))
	zostogafile = os.path.join(outdir, "{}_ZOSTOGA.pkl".format(args.pipeline_id))
	zosfile = os.path.join(outdir, "{}_ZOS.pkl".format(args.pipeline_id))
	tlmfile = os.path.join(outdir, "{}_tlmdata.pkl".format(args.pipeline_id))

	# Run the OD preprocessing if intermediate files are not present
	scenario_dsl = args.scenario_dsl
	if len(scenario_dsl) == 0:
		scenario_dsl = args.scenario

	if os.path.isfile(configfile) and os.path.isfile(zostogafile) and os.path.isfile(zosfile):
		print("{}, {}, and {} found, skipping OD preprocessing".format(configfile, zostogafile, zosfile))
	else:
		tlm_preprocess_oceandynamics(scenario_dsl, args.model_dir, not args.no_drift_corr, args.no_correlation, args.pyear_start, args.pyear_end, args.pyear_step, args.locationfile, args.baseyear, args.pipeline_id)

	# Done
	sys.exit()
