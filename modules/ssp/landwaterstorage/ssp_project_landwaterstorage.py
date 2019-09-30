import numpy as np
#import ntpath
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.special import erf
import argparse
import os
import pickle
import time
from netCDF4 import Dataset

''' ssp_project_landwaterstorage.py

Code generated 16-09-2019, by Tim Hermans

This script runs the land water storage projections for the SSP LWS workflow. This task 
projects the future contribution of groundwater depletion and reservoir impoundment
to global mean sea level based on the selected SSP of population growth.

Parameters:
Nsamps = Number of samples to project
rng_seed = Seed value for the random number generator
pipeline_id = Unique identifier for the pipeline running this code

Output:
"%PIPELINE_ID%_projections.pkl" = Pickle file that contains the global LWS projections

'''


def ssp_project_landwaterstorage(Nsamps, rng_seed, pipeline_id):
	
	# Load the fit file
	fitfile = "{}_fit.pkl".format(pipeline_id)
	try:
		f = open(fitfile, 'rb')
	except:
		print("Cannot open fit file\n")
	
	# Extract the fit variables
	my_fit = pickle.load(f)
	f.close()
	
	popscen = my_fit['popscen']
	popscenyr = my_fit['popscenyr']
	#popscenids = my_fit['popscenids']
	dams_popt = my_fit['dams_popt']
	mean_dgwd_dt_dpop = my_fit['mean_dgwd_dt_dpop']
	std_dgwd_dt_dpop = my_fit['std_dgwd_dt_dpop']
	
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
	dam_pcterr = my_config['dam_pcterr']
	yrs = my_config['yrs']
	scen = my_config['scen']
	dotriangular = my_config['dotriangular']
	includepokhrel = my_config['includepokhrel']
	
	# optimisation problem, least squares of fitting dams with sigmoidal function of population
	def sigmoidal(pop0,a,b,c,I0):
		return(a*erf((pop0/1e6-b)/c) + I0) # see Kopp et al. 2014 eq.1

	##################################################
	#select scenario population using target RCP or SSP scenario	
	# prefered SSP RCP combinations (correspondence with Aimee)    
	RCPtoSSP = {
		"rcp19": "ssp1",
    	"rcp26": "ssp1",
        "rcp45": "ssp2",
        "rcp70": "ssp3",
        "rcp85": "ssp5",
        }     
        
    # SSP ordered from low to high projections
	SSPorder = {
        "ssp1": 0,
        "ssp5": 1,
        "ssp2": 2,
        "ssp4": 3,
        "ssp3": 4,
        }  
    
    # extract SSP scenario from configured target RCP or SSP scenario
	if scen[0:3] == "rcp":
		targetSSP = RCPtoSSP[scen]
		if scen not in RCPtoSSP:
			raise Exception('Configured RCP scenario does not have a preferred SSP combination.')
	else:
		targetSSP = scen        
    
	# draw scenario population from target scenario 
	popdraw = popscen[:,SSPorder[targetSSP]]
	
	# interpolate to annual means
	popdraw = np.interp(np.linspace(2000,2100,101),popscenyr,popdraw)
	popscenyr = np.linspace(2000,2100,101)

	# random draw functions for reservoirs and gwd

	# symbolic function interpolating the cumulative sum of population for a
	# certain probability x1 times the normal distribution with mean dgwd/dt/dpopmGWD
	# and standard deviation stdGWD defined earlier, onto desired years    
	if dotriangular == 0:
		def gwddraw(seed1):
			return( np.interp( yrs, popscenyr, np.cumsum(np.multiply(popdraw,\
				(mean_dgwd_dt_dpop + norm.ppf(seed1)*std_dgwd_dt_dpop))) ) )
	else:
		def gwddraw(seed1,seed2):
			return( np.interp( yrs, popscenyr, np.cumsum(np.multiply(popdraw,\
				np.interp(seed1, np.array([0, .5, 1]), np.array([std_dgwd_dt_dpop[0], mean_dgwd_dt_dpop,std_dgwd_dt_dpop[-1]]))))\
				* (1 + norm.ppf(seed2)*dgwd_dt_dpop_pcterr) ) )

	# random draw from reservoir storage: sigmoidal function with 
	# randomly drawn population as input (>t=2000, see Kopp 2014)
	# and subtracts the sigmoidal function with the population at t=2000 (this
	# will then be the origin). 

	# This is multiplied with a normal distribution with a mean of 1 and a std of 
	# 25% (default defined error elated to the impoundment rate. Kopp 2014: 2sigma=50%):
	# - minus sign since reservoir storage leads to GSL drop - 

	pop2000 = pop0[t0==2000]  #population at 2000

	def damdraw(seed1): ###decide max of sigmoid or pop2000->choose
		poprand = np.array([popdraw]) # random draw from population
		poprand[poprand<pop2000] = pop2000 # impoundment is not allowed to be reduced below yr 2000 levels (Kopp et al., 2014)

		X = np.multiply(-1*(sigmoidal(poprand,dams_popt[0],dams_popt[1],dams_popt[2],dams_popt[3])\
							- sigmoidal(pop2000[0],dams_popt[0],dams_popt[1],dams_popt[2],dams_popt[3])), 1+norm.ppf(seed1)*dam_pcterr )
    
		#X = np.multiply(-1*(sigmoidal(np.array([popdraw(seed1)]),dams_popt[0],dams_popt[1],dams_popt[2],dams_popt[3])\
		#                    - sigmoidal(pop2000[0],dams_popt[0],dams_popt[1],dams_popt[2],dams_popt[3])), 1+norm.ppf(seed2)*dam_pcterr )
    
		return( np.interp(yrs, popscenyr, X[0] ) )

	##################################################
	# generate seeds and draw samples 
	np.random.seed(rng_seed)
	if isinstance(Nsamps,int): #if only given a single number Nsamps
		seeds0 = np.linspace(0,1,Nsamps+2)
		seeds0 = seeds0[1:-1]
	else:
		seeds=Nsamps
		seeds0=Nsamps

	if seeds0.ndim == 1: # if seeds is a vector
		seeds = np.empty((4,len(seeds0)))
		for j in range(0,4):
			seeds[j,:] = seeds0[np.random.permutation(len(seeds0))]
        
	#draw the samples
	damsamps = np.empty((len(yrs),Nsamps))   
	gwdsamps = np.empty((len(yrs),Nsamps)) 

	if dotriangular == 0:       
		for ii in range(0,np.shape(seeds)[1]):
			damsamps[:,ii] = damdraw(seeds[3,ii])
			gwdsamps[:,ii] = gwddraw(seeds[1,ii])
	else:
		for ii in range(0,np.shape(seeds)[1]):
			damsamps[:,ii] = damdraw(seeds[3,ii])
			gwdsamps[:,ii] = gwddraw(seeds[1,ii],seeds[2,ii])
    
	# add to total lws equivalent gsl    
	lwssamps = gwdsamps + damsamps
	
	# Store the variables in a pickle
	output = {'lwssamps': lwssamps, 'years': yrs, 'scen': scen}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()
	
	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(yrs))
	samp_dim = rootgrp.createDimension("samples", Nsamps)

	# Populate dimension variables
	year_var = rootgrp.createVariable("year", "i4", ("years",))
	samp_var = rootgrp.createVariable("sample", "i8", ("samples",))

	# Create a data variable
	samps = rootgrp.createVariable("samps", "f4", ("years", "samples"), zlib=True, least_significant_digit=2)
	
	# Assign attributes
	rootgrp.description = "Global SLR contribution from land water storage according to Kopp 2014 workflow"
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0} - {1}".format(pipeline_id, scen)
	year_var.units = "[-]"
	samp_var.units = "[-]"
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = yrs
	samp_var[:] = np.arange(0,Nsamps)
	samps[:,:] = lwssamps

	# Close the netcdf
	rootgrp.close()
	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the land water storage projection stage for the SSP LWS SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the SSP LWS module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', '-n', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', '-s', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the preprocessing stage with the provided arguments
	ssp_project_landwaterstorage(args.nsamps, args.seed, args.pipeline_id)
	
	exit()