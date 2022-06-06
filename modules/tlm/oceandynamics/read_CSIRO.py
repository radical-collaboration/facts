import numpy as np

''' read_CSIRO.py

This script loads the global sea level data from CSIRO file

Parameters: 
file = Full path to file for loading
centeryear = Year at which the data should be centered (default: 2005)

Return: 
X1 = Matrix containing latitude, longitude, and years of record (lat and lon set to 1e6)
Y = Sea level data
dY = Uncertainty around sea level data
regions = Region identifer
regionsu = ?
sitenames = "CSIRO GSL"
sitecoords = [1e6 1e6]
sitelen = Length of the global sea level record
sitecoastline = 0

Note: This function returns many fields that are not necessary when dealing with a global
data set. The reason for having them though is to maintain commonalities with the full
ReadPSMSLData.m code from the original K14 workflow.

'''

def read_CSIRO(file, centeryear=2005):
	
	# Open the file
	data = np.loadtxt(file, skiprows=1, delimiter=',')
	
	# Extract the data
	years = data[:,0]
	rsl = data[:,1]
	rslunc = data[:,2]/2
	
	# Center the data on the required year
	centeryear_ind = np.argmin(np.abs(years - centeryear))
	rsl = rsl - rsl[centeryear_ind]
	
	# Populate the output variables
	X1 = np.hstack((1e6*np.ones((len(years),2)), years.reshape((-1,1))))
	regions = np.zeros(len(years))
	Y = rsl
	dY = rslunc
	sitelen = len(rsl)
	regionsu = 0
	sitenames = "CSIRO GSL"
	sitecoords = np.array([1e6,1e6])
	sitecoastline = 0
	
	return(X1,Y, dY, regions, regionsu, sitenames, sitecoords, sitelen, sitecoastline)