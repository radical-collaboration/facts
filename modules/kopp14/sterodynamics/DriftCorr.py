import numpy as np
import pickle
import sys
from read_CSIRO import *
from Smooth import *
from readMarzeion import *

''' DriftCorr.py

This script applies a drift correction to ZOSTOGA

Parameters: 
ZOSTOGA = Time series of thermosteric sea-level change for each included model as derived
		  by the IncludeModels.py or the SmoothZOSTOGA.py script (years)
years = Years of interest
baseyear = Year that is used to center the ZOSTOGA and global sea-level values
scenario = RCP scenario (e.g "rcp85")
gslfile = Full path to Global sea-level rise data


Return: 
dcZOSTOGA = Drift-corrected global average thermosteric sea-level change (years)
CWdrift = Drift in global sea-level change
histGICrate = Historical GIC rate
selectyears = Years that coincide with years and the GIC records

'''

def DriftCorr(ZOSTOGA, years, baseyear, scenario, gslfile):
	
	# Read in the global SLR data
	(GSLx,GSLy,_,_,_,_,_,_,_) = read_CSIRO(gslfile)
	
	# Smooth the sea-level signal and convert to meters
	GSLy = Smooth(GSLy, 19)/1000
	
	# Calculate the CWDrift value
	inds = np.nonzero(GSLx[:,2] <= 1900)[0]
	CWdrift = (GSLy[inds.max()] - GSLy[inds[0]]) / (GSLx[inds.max(),2]-GSLx[inds[0],2])
	
	# Center the SLR on baseyear
	baseyear_ind = np.nonzero(np.floor(GSLx[:,2])==baseyear)[0]
	GSLy = GSLy - np.mean(GSLy[baseyear_ind])
	
	# Determine the years over which the corrections are applied
	selectyears = np.array([1861, np.floor(GSLx[inds.max(),2])])
	
	# Load the glacier data
	glacdir = os.path.join(os.path.dirname(__file__), "Marzeion2012supplement")
	fpmap = os.path.join(os.path.dirname(__file__), "fingerprint_region_map.csv")
	(projGIC85, projGIC85se, projGIC85yrs, projGIC85model,_,_,_) = readMarzeion(scenario, glacdir, fpmap, discardAntarctica=True)
	
	# Find which indices in the glacier data correspond to the start and end years in selectyears
	startyear_ind = np.nonzero(projGIC85yrs[:,0] == selectyears[0])[0]
	endyear_ind = np.nonzero(projGIC85yrs[:,0] == selectyears[1])[0]
	
	# Sum the contributions across regions
	histGIC = np.sum(projGIC85, axis=1)
	
	# Calculate the historic rate of change of GIC contributions and convert to meters
	histGICrate = (histGIC[endyear_ind,:] - histGIC[startyear_ind,:]) / np.diff(selectyears) / 1000
	
	# Find the difference between the CWdrift and the mean historical GIC rate
	# Note: Units are in meters
	histresrate = CWdrift - np.mean(histGICrate)
	
	# Find where in the ZOSTOGA data the years that match selectyears
	ZOSTOGA_startyear_ind = np.nonzero(years == selectyears[0])[0]
	ZOSTOGA_endyear_ind = np.nonzero(years == selectyears[1])[0]
	
	# Calculate the ZOSTOGA drift
	drift = (ZOSTOGA[ZOSTOGA_endyear_ind,:] - ZOSTOGA[ZOSTOGA_startyear_ind,:]) / np.diff(selectyears)
	
	# Apply the pure ZOSTOGA drift correction
	ZOSTOGA_tmp = ZOSTOGA - np.outer(years - selectyears[0], drift)
	
	# Apply the drift correction from the GSL and glacier data
	dcZOSTOGA = ZOSTOGA_tmp + (histresrate * (years - selectyears[0]))[:,np.newaxis]
	
	# Center the drift-corrected ZOSTOGA to the baseyear
	ZOSTOGA_baseyear_ind = np.nonzero(years == baseyear)[0]
	dcZOSTOGA = dcZOSTOGA - dcZOSTOGA[ZOSTOGA_baseyear_ind,:]
	
	# Return variables
	# Note: dcZOSTOGA, CWdrift, and histGICrate are all in meters
	return(dcZOSTOGA, CWdrift, histGICrate, selectyears)

if __name__ == '__main__':
	
	gslfile = os.path.join(os.path.dirname(__file__), "CSIRO_Recons_gmsl_yr_2011.csv")
	
	# Load the ZOSTOGA file
	zostogafile = os.path.join(os.path.dirname(__file__), "kopp14_thermalexp_ZOSTOGA.pkl")
	try:
		f = open(zostogafile, 'rb')
	except:
		print("Cannot open ZOSTOGA file\n")
	
	# Extract the configuration variables
	my_zostoga = pickle.load(f)
	f.close()

	ZOSTOGA = my_zostoga["ZOSTOGA"]
	
	x = DriftCorr(ZOSTOGA,np.arange(1861,2300),2000, "rcp85", gslfile)
	#print(x)
	
	exit()