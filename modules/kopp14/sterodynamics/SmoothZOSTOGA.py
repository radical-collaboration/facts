import numpy as np
import sys
from Smooth import Smooth

''' SmoothZOSTOGA.py

This script fixes the ZOSTOGA "suturing" problem and applies smoothing.

Parameters: 
ZOSTOGA = Time series of thermosteric sea-level change for each included model as derived
		  by the IncludeModels.py script (years, nmodels)
years = Years of interest
baseyear = Base year from which to center ZOSTOGA on the mean
smoothwin = Number of years contained in the smoothing window

Return: 
ZOSTOGA = Global average thermosteric sea-level change with suturing (years, nmodels)
sZOSTOGA = Smoothed ZOSTOGA (years, nmodels)

'''

def SmoothZOSTOGA(ZOSTOGA, years, baseyear, smoothwin):
	
	# Initialize the smoothed ZOSTOGA variable
	sZOSTOGA = np.nan * ZOSTOGA
	
	## Center ZOSTOGA on its mean
	# Find the indices that are within 10 years of the baseyear
	temp1 = years <= (baseyear+10)
	temp2 = years >= (baseyear-10)
	center_inds = np.flatnonzero(temp1 * temp2)
	good_inds = center_inds[np.flatnonzero(~np.isnan(ZOSTOGA[center_inds]))]
	if(len(good_inds) > 0):
		ZOSTOGA = ZOSTOGA - np.mean(ZOSTOGA[good_inds])
	
	# Patch the suturing problem
	syear_ind = int(np.flatnonzero(years == 2007))
	diffa = (ZOSTOGA[syear_ind] - ZOSTOGA[syear_ind-1])
	diffb = (ZOSTOGA[syear_ind-2] - ZOSTOGA[syear_ind-1])
	if(np.abs(diffa) > 20*np.abs(diffb)):
		offset = diffa - diffb
		ZOSTOGA[syear_ind:] = ZOSTOGA[syear_ind:] - offset
	
	# Smooth ZOSTOGA
	good_inds = np.nonzero(~np.isnan(ZOSTOGA))[0]
	sZOSTOGA[good_inds] = Smooth(ZOSTOGA[good_inds], smoothwin)
	
	# Center the smoothed ZOSTOGA on the baseyear
	sZOSTOGA = sZOSTOGA - sZOSTOGA[np.flatnonzero(years == baseyear)]
	
	return(ZOSTOGA, sZOSTOGA)