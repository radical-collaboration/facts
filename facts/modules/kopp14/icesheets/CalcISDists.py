import numpy as np
from FitLNDistQuants import FitLNDistQuants as fitLN

'''
CalcISDists

Fit log-normal distributions to the Bamber and Aspinall (2013) and IPCC AR5 ice sheet
melt estimates.

Parameters:
barates2100 = BA13 ice sheet melt rates in 2100
lastdecadegt = Ice sheet mass loss in giga-tonnes
aris2090 = IPCC AR5 ice sheet contributions to SLR in 2090

Return:
batheteais = Distribution parameters for east-antarctic ice sheet from BA13 rates
bathetwais = Distribution parameters for west-antarctic ice sheet from BA13 rates
bathetgis = Distribution parameters for Greenland ice sheet from BA13 rates
arthetais = Distribution parameters for antarctic ice sheet from AR5 estimates
arthetgis = Distribution parameters for the Greenland ice sheet from AR5 estimates
islastdecade = Ice sheet contributions to SLR over last decade in SLR-equivalent

Note: This function is transcoded from Robert Kopp's CalculateISDists.m.
'''

def CalcISDists(barates2100, lastdecadegt, aris2090):
	
	# Use Shepherd estimates for last decade (2000 - 2011), pooling APIS and WAIS
	surfacearea = 3.16E14  # Surface area of ice sheet in m^2
	densitysw = 1020  # Density of (salt?) water in kg/m^3
	gt2mm = -1E12/surfacearea/densitysw * 1000  # Convert giga-tonnes to mm
	
	# Convert the mass loss of ice sheets over last decade from giga-tonnes to mm
	# equivalent of sea-level rise
	islastdecade = lastdecadegt * gt2mm
	
	# Fit log-normal distributions to the BA13 rates
	batheteais = fitLN(barates2100[2,1], barates2100[2,2:5], barates2100[2,0], [0.05,0.5,0.95])
	bathetwais = fitLN(barates2100[1,1], barates2100[1,2:5], barates2100[1,0], [0.05,0.5,0.95])
	bathetgis = fitLN(barates2100[0,1], barates2100[0,2:5], barates2100[0,0], [0.05,0.5,0.95])
	
	# Combine the west and east antarctic ice sheet contributions over last decade
	arislastdecade = np.array([islastdecade[0], islastdecade[1]+islastdecade[2]])
	
	# Calculate a baseline rate for the AR5 estimate
	arbase = (2090 - 1995) * arislastdecade
	
	# Normalize the ice sheet rate estimates from AR5 to a common baseline so that the
	# rate for year 2100 can be estimated
	arisaddl0 =  aris2090 - arbase.reshape(2,1)
	arisaccel0 = 2*arisaddl0/(2090-2011)**2
	
	# Setup the rates matrix for the AR5 estimates
	arrates2100 = (arisaccel0*(2100-2011)) + arislastdecade.reshape(2,1)
	
	# Fit log-normal distributions to the AR5 estimates
	arthetais = fitLN(arrates2100[1,1], arrates2100[1,:], -3, [0.167,0.5,0.833])
	#arthetais = np.array([-169.3219, 5.1349, 0.0156])   # Matlab version's results
	arthetgis = fitLN(arrates2100[0,1], arrates2100[0,:], 0, [0.167,0.5,0.833])
	
	return(batheteais, bathetwais, bathetgis, arthetais, arthetgis, islastdecade)