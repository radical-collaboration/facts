import numpy as np
from numpy import matlib as mb
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.special import erf
import argparse
import pickle
import os

''' ssp_fit_landwaterstorage.py

Code generated 16-09-2019, by Tim Hermans

This script runs the land water storage fitting task for the SSP LWS workflow. 
This task fits a linear relation between historical population data and groundwater 
depletion, and a sigmoidal relation between population and reservoir impoundment. 

Parameters: 
pipeline_id = Unique identifier for the pipeline running this code

Output:
"%PIPELINE_ID%_fit.pkl" = Pickle file that contains the fitted submodel information

'''

def extend_pop(popscen, popscenyrs):
	
	# Rate as a function of population scenario (in percent)
	# Obtained from https://www.un.org/development/desa/pd/sites/www.un.org.development.desa.pd/files/files/documents/2020/Jan/un_2002_world_population_to_2300.pdf
	# Table 1, page 14 (pdf page 28)
	rates = np.array([[0.76, 0.04, -0.46, -0.74, -0.75, -0.60, -0.48, -0.38, -0.32, -0.31, -0.31, -0.32],
						[1.03, 0.51, 0.13, -0.07, -0.15, -0.11, -0.03, 0.03, 0.06, 0.06, 0.05, 0.05], 
						[1.28, 0.96, 0.64, 0.46, 0.35, 0.36, 0.45, 0.51, 0.54, 0.54, 0.54, 0.54]]) / 100 + 1
	
	# Extend these rates over the 25-year period they represent
	ext_rates = np.repeat(rates, 25, axis=1)
	ext_years = np.arange(2001, 2301)	# Rates are for years 2001 - 2300
	
	# SSP map to population scenario
	# "popscen" is sorted in terms of lowest to highest projected population, which is
	# SSP1, SSP5, SSP2, SSP4, SSP3
	scenario_map = np.array([0,0,1,1,2])
	
	# Select the appropriate rates
	scenario_rates = ext_rates[scenario_map,:].T
	
	# Which year in the extended years matches the next population year
	# Note: SSP projection data should end in year 2100
	year_idx = np.flatnonzero(ext_years == popscenyrs[-1]) + 1
	
	# Initialize variables to hold the extended population and years
	# These will be appended to popscen and popscenyrs provided
	ext_pop = []
	ext_yrs = []
	
	# Set the initial population
	pop0 = popscen[-1,:]
	
	# Loop over the extended years
	for i in np.arange(year_idx, len(ext_years)):
		
		# Calculate the population for this year
		ext_pop.append(pop0 * scenario_rates[i,:])
		
		# Set the current population for the next iteration
		pop0 = np.array(ext_pop[int(i - year_idx)])
		
		# Append this year to the extended years
		ext_yrs.append(ext_years[i])
	
	# Cast boundary rates and years from a list into a numpy arrays
	ext_pop = np.array(ext_pop)
	ext_yrs = np.array(ext_yrs)
	
	# Append the extension onto the original data
	popscen = np.concatenate((popscen, ext_pop), axis=0)
	popscenyrs = np.concatenate((popscenyrs, ext_yrs))
	
	return(popscen, popscenyrs)	


def ssp_fit_landwaterstorage(pipeline_id):
	
	# Load the data file
	datafile = "{}_data.pkl".format(pipeline_id)
	try:
		f = open(datafile, 'rb')
	except:
		print("Cannot open data file\n")
	
	# Extract the data variables
	my_data = pickle.load(f)
	f.close()
	
	t = my_data["t"]
	pop = my_data["pop"]
	tdams = my_data["tdams"]
	dams = my_data["dams"]
	tgwds = my_data["tgwd"]
	gwds = my_data["gwd"]
	popscenyr = my_data["popscenyr"]
	popscen = my_data["popscen"]
	
	# Load the configuration file
	configfile = "{}_config.pkl".format(pipeline_id)
	try:
		f = open(configfile, 'rb')
	except:
		print("Cannot open config file\n")
	
	# Extract the configuration variables
	my_config = pickle.load(f)
	f.close()
	
	t0 = my_config['t0']
	pop0 = my_config['pop0']
	dgwd_dt_dpop_pcterr = my_config['dgwd_dt_dpop_pcterr']
	yrs = my_config['yrs']
	scen = my_config['scen']
	dotriangular = my_config['dotriangular']
	includepokhrel = my_config['includepokhrel']



	# interpolate reservoir years to population history years t0
	dams = np.interp(t0,tdams,dams)

	# optimisation problem, least squares of fitting dams with sigmoidal function of population
	def sigmoidal(pop0,a,b,c,I0):
		return(a*erf((pop0/1e6-b)/c) + I0) # see Kopp et al. 2014 eq.1

	# initial guess
	pinit = np.array([ max(dams)/2, 1, 1, max(dams)/2 ])
	# curve fit
	dams_popt, dams_pcov = curve_fit(sigmoidal, pop0, dams,p0=pinit)
	
	# Initialize variables
	dgwd_dt_dpop_all = [] # change of gwd rate w.r.t. population
	dgwd_dt_all      = [] # gwd rate
	pop2gwd_all      = [] # population
	
	# Loop over each of the GWD files' data
	for i in np.arange(0,2+includepokhrel):
		
		# Working with these particular data
		gwd = gwds[i,:]
		tgwd = tgwds[i,:]
		
		# Remove entries where either tgwd or gwd are "nan"
		inds = np.intersect1d(np.flatnonzero(~np.isnan(tgwd)), np.flatnonzero(~np.isnan(gwd)))
		tgwd = tgwd[inds]
		gwd = gwd[inds]

		# remove duplicate years        
		tgwd, ui = np.unique(tgwd, return_index=True)
		gwd = gwd[ui]
		
		# interpolate to population history years t (5 year samples), without etrapolating outside bounds
		bound = np.where((tgwd[0] <= t) & (t <= tgwd[-1]))
		bound = bound[0]
		gwd = np.interp(t[bound],tgwd,gwd)
		
		dgwd_dt = np.diff(gwd)/np.diff(t[bound]) # compute dgwd/dt
		pop2gwd = (pop[bound[0:len(bound)-1]]+pop[bound[1:len(bound)]])/2 # get the population at the timesteps between timesteps at which there is GWD information
		dgwd_dt_dpop = np.linalg.lstsq(np.transpose([pop2gwd]),np.transpose([dgwd_dt]),rcond=None)  # solve least squares linear fit for the change of ground depletion rate with population
		
		# store dgwd_dt, dgwd_dt_dpop, pop2gwd
		dgwd_dt_dpop_all = np.concatenate((dgwd_dt_dpop_all,dgwd_dt_dpop[0][0]))
		dgwd_dt_all.append(dgwd_dt)
		pop2gwd_all.append(pop2gwd)
    
	# mean and standard deviation dgwd/dt/dpop
	if dotriangular == 0: # if no triangular distribution
		mean_dgwd_dt_dpop = np.mean(dgwd_dt_dpop_all)
		std_dgwd_dt_dpop = np.sqrt( np.std((dgwd_dt_dpop_all),ddof=1)**2 + (mean_dgwd_dt_dpop*dgwd_dt_dpop_pcterr)**2)
	
	else: # if triangular distribution
		mean_dgwd_dt_dpop = np.median(dgwd_dt_dpop_all)
		std_dgwd_dt_dpop = np.array((min(dgwd_dt_dpop_all),max(dgwd_dt_dpop_all)))


	popscenyr0 = popscenyr
	popscen0 = popscen
	popscen = np.divide(popscen,1e3) # convert to thousands
	#popscen = popscen/1000.0

	# if scenarios start >2000, prepend the years from 2000 onward and
	# interpolate historical population up to the start of projections
	if min(popscenyr)>2000:
		popscenyr = np.insert(popscenyr,0,range(2000,int(min(popscenyr))))
		popscen = np.concatenate((np.transpose(mb.repmat(np.interp(range(2000,int(min(popscenyr0))),t0,pop0),5,1)),popscen))
	
	# Extend the population projections
	(popscen, popscenyr) = extend_pop(popscen, popscenyr)
	
	# Test to make sure the target years are within the projected population years
	if max(yrs)>max(popscenyr):
		raise Exception("Target year {} exceeds the maximum year of SSP population scenarios {}".format(max(yrs), max(popscenyr)))
	
	# Store the variables in a pickle
	output = {'popscen': popscen, 'popscenyr': popscenyr,\
				'dams_popt': dams_popt, 'mean_dgwd_dt_dpop': mean_dgwd_dt_dpop, \
				'std_dgwd_dt_dpop': std_dgwd_dt_dpop}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_fit.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()


if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the land water storage fitting stage for the SSP LWS SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the SSP LWS module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ssp_fit_landwaterstorage(args.pipeline_id)
	exit()