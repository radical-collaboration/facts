import numpy as np

''' DriftCorr.py

This script applies a drift correction to ZOSTOGA

Parameters: 
ZOSTOGA = Time series of thermosteric sea-level change for each included model as derived
		  by the IncludeModels.py or the SmoothZOSTOGA.py script (years)
years = Years of interest
baseyear = Year that is used to center the ZOSTOGA and global sea-level values
modellist = List of included models produced by IncludeModels.py


Return: 
dcZOSTOGA = Drift-corrected global average thermosteric sea-level change (years)
CWdrift = Drift in global sea-level change
histGICrate = Historical GIC rate
selectyears = Years that coincide with years and the GIC records

'''

def DriftCorr(ZOSTOGA, years, baseyear, modellist):
	
	# Initialize variables
	dcZOSTOGA = ZOSTOGA
	CWdrift = np.nan
	histGICrate = np.nan
	selectyears = np.nan
	
	return(dcZOSTOGA, CWdrift, histGICrate, selectyears)