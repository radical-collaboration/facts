import os
import numpy as np
import pandas as pd
import pickle
import argparse
from scipy.stats import genpareto
from gplike import gplike

# Do not warn about chained assignments
pd.options.mode.chained_assignment = None  # default='warn'

""" extremesl_fit.py
This script fits a General Pareto Distribution on preprocessed GESLA2 data. It
also calculates the associated covariance matrix.

The preprocessing is based on the MATLAB code of Thomas Frederikse used for SROCC.

Input parameters:
station_data        = dictionary with relevant preprocessed GESLA 2 data
gdp_pot_treshold    = user configured threshold

Output: adds GDP parameters to dictionary with station data
loc                 = location parameter (POT threshold)
gp_shape            = shape factor
gp_scale            = scale factor
gp_cov              = GPD parameter covariances
extremes_loc        = exremes exceeding the location parameter
avg_exceed          = average exceedances of threshold per year
decl_esls_pyear     = number of declustered extremes in each year
esl_years           = years with extremes

Created on Mon Oct 30 16:43:56 2019
@author: Tim Hermans
tim(dot)hermans@nioz(dot)nl

Updated on Tue Apr 28 2020 at 15:13 EDT
Adapted for use in FACTS by Gregory Garner
gregory(dot)garner@rutgers(dot)edu

"""

def extremesl_fit(station_data_file, pipeline_id):
	
	# Load the station data file
	try:
		f = open(station_data_file, 'rb')
	except:
		print("Cannot open station data file: {}\n".format(station_data_file))
	
	# Extract the configuration variables
	station_data = pickle.load(f)
	f.close()
	
	# Keep track of station IDs that don't produce valid GPD fit
	bad_stations = []
	
	# Initialize dictionary to hold fitted data
	fitted_data = {}
	
	# Loop through the available station data
	for station_id in station_data:
		
		#range of Peak-Over-Threshold values and corresponding heights at station
		potvals = station_data[station_id]['pot_vals']
		
		#etreme values based on POT with 95%
		extremes_prepo = station_data[station_id]['decl_extremes']
	
		#Determine location parameter, i.e. find the peak height corresponding to POT threshold in projection settings
		#location parameter = threshold water-level above which return levels are estimated with the GPD,
		#loc = potvals[np.where(potvals[:,0] == gpd_pot_threshold)][0][1]
		loc = potvals
		
		#extremes that exceed location parameter
		extremes_loc = extremes_prepo[extremes_prepo['height'] > loc]
		extremes_loc_timesorted = extremes_loc.sort_values(by=['mytime'])
	
		#count declustered days with extreme obs per year
		decl_esls_pyear = extremes_loc.groupby(extremes_loc_timesorted.mytime.dt.year,as_index=True).count()
		esl_years = np.unique(extremes_loc_timesorted['mytime'].dt.year.values)
		
		#average exceedances per year for configured threshold
		avg_exceed = 365.25*24*len(extremes_loc)/len(station_data[station_id]['obs'])
		#^ above is not completely correct, because there are less than so many hours in one year? should take the average number of hours of obs in one year?
	
		##### FIT GPD
		#provide initial guess using method of moments (code based on gpfit.m)
		xbar = np.mean(extremes_loc['height'].values-loc) #mean of extremes
		s2 = np.var(extremes_loc['height'].values-loc) #variance of extremes
		k0 = -.5 * ((xbar**2)/s2 - 1) #initial guesses
		sigma0 = .5 * xbar * ( (xbar**2) / s2 + 1)
		xmax = max(extremes_loc['height'].values-loc)
		
		# Method of moments invalid (code based on gpfit.m)
		if (k0 < 0 and xmax >= -sigma0/k0 ):
			#assume exponential distribution
			k0 = 0
			sigma0 = xbar
		
		#fit gpd based on exceedences and initial guess, note that optimization vals differ slightly from gpfit in MATLAB
		gp_params = genpareto.fit(extremes_loc['height'].values-loc,loc=0,scale=sigma0)
		
		#calculate covariance matrix of estimated parameters
		gp_nlogl, gp_cov = gplike(gp_params[0], gp_params[2],extremes_loc['height'].values-loc)
		
		# Are the GPD parameters outside the supported range? If so, add this station to the list of bad stations.
		if gp_nlogl == np.Inf:
			print('GDP parameters of '+station_id+' are outside of the supported range, confidence intervals and standard errors cannot be computed reliably.')
			print(station_id+' will be ommitted from the results.')
			bad_stations.append(station_id)
		
		###### QUERY HISTORICAL RETURN FREQUENCIES FOR TEST HEIGHTS
		#get MHHW for Gumbel distribution below Pareto
		obs=station_data[station_id]['obs']
		obs2D = obs.groupby(pd.Grouper(key='mytime',freq='2D')).max() #2-daily maxima
		mhhw = np.nanmean(obs2D.height) #calculate MHHW as mean of maximum value per 2 days (to account for tidal cycle longer than 1 day)
		mhhwFreq = 365.25/2; #assume MHHW is exceeded once every 2 days
		
		##### STORE DATA IN DICTIONARY
		fitted_data[station_id] = {'loc': loc, 'gp_shape': gp_params[0], 'gp_scale': gp_params[2], \
			'gp_cov': gp_cov, 'extremes_loc': extremes_loc, 'avg_exceed': avg_exceed, 'decl_esls_pyear': decl_esls_pyear, \
			'esl_years': esl_years, 'mhhw': mhhw, 'mhhwFreq': mhhwFreq, 'obs': obs}
			
		#fitted_data[station_id]['gp_shape'] = gp_params[0]
		#fitted_data[station_id]['gp_scale'] = gp_params[2]
		#fitted_data[station_id]['gp_cov'] = gp_cov
		#fitted_data[station_id]['extremes_loc'] = extremes_loc
		#fitted_data[station_id]['avg_exceed'] =avg_exceed
		#fitted_data[station_id]['decl_esls_pyear'] =decl_esls_pyear
		#fitted_data[station_id]['esl_years'] = esl_years
		#fitted_data[station_id]['mhhw'] = mhhw
		#fitted_data[station_id]['mhhwFreq'] = mhhwFreq
	
	# Remove the stations that don't support GPD fit
	for bad_station in bad_stations:
		del fitted_data[bad_station]
	
	# If there are no more stations left, raise an error 
	if len(fitted_data) == 0:
		raise Exception("No stations available with valid GPD fit")
	
	# Write station data to an output pickle
	outdir = os.path.dirname(__file__)
	outfile = open(os.path.join(outdir, "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(fitted_data, outfile)
	outfile.close()
	
	return(0)



if __name__ == "__main__":
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the fitting stage for the extreme sea-level workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--station_data_file', help="Preprocessed station data from GESLA database (produced from preprocessing stage)", default=None)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Use default for station data file if necessary
	if args.station_data_file is None:
		station_data_file = os.path.join(os.path.dirname(__file__), "{}_station_data.pkl".format(args.pipeline_id))
	else:
		station_data_file = args.station_data_file
	
	# Run the fitting process on the 
	extremesl_fit(station_data_file, args.pipeline_id)
	
	exit()