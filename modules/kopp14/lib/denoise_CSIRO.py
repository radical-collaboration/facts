import numpy as np
import numpy.linalg as lin
import os
from covarfuns import *

''' denoise_CSIRO.py

This script denoises the CSIRO data or data read in from the "read_CSIRO.py" script

Parameters: 
thetTG = Theta values for the target sites (in this context, global)
thinyrs = Thin the data set to years at this interval
firstyear = First year desired in the data set

Return: 


Note: 

'''

def denoise_CSIRO(TGcoords, TGrsl, TGrslunc, thetTG, thinyrs=1, firstyear=np.nan):
	
	# Define covariance function for TG
	def cvfuncTG(t1,t2,dt1t2,thetas,ad):
		term1 = kDP(t1,t2,1) * (thetas[0]**2 + kDELTAG(ad,thetas[6]) + kGEOG(ad,thetas[9:11]))
		term2 = kRQ(dt1t2, np.append(1, thetas[2:4])) * (thetas[1]**2 + kDELTAG(ad,thetas[7]) + kGEOG(ad,thetas[11:13]))
		term3 = kMat1(dt1t2,np.append(1,thetas[5])) * (thetas[4]**2 + kDELTAG(ad,thetas[8]) + kGEOG(ad, thetas[13:15]))
		return(term1 + term2 + term3)
	
	# Get the years of the data from TGcoords...
	TGyears = TGcoords[:,2]
	
	# Determine firstyear if necessary
	if(np.isnan(firstyear)):
		firstyear = TGyears[0]
	
	# Get the target years.  Apply thinning if necessary
	TGtargyears = np.arange(firstyear, TGyears[len(TGyears)-1]+1, thinyrs)
	
	# Generate the target coordinates of same length as the target years
	TGtargcoords = np.full((len(TGtargyears), 2), 1e6)
	
	# Append the coordinates and years together...
	TGtargcoords = np.hstack((TGtargcoords, TGtargyears.reshape(-1,1)))
	
	# Define training and testing covariance functions for GPR
	def traincv(thet):
		return(cvfuncTG(TGyears,TGyears,dYears(TGyears,TGyears),thet,dDist(TGcoords,TGcoords)) + np.diag(TGrslunc**2))
	def testcv(thet):
		return(cvfuncTG(TGyears,TGtargyears,dYears(TGyears,TGtargyears),thet,dDist(TGcoords,TGtargcoords)))
	def testcv2(thet):
		return(cvfuncTG(TGtargyears,TGtargyears,dYears(TGtargyears,TGtargyears),thet,dDist(TGtargcoords,TGtargcoords)))
		
	# Noise mask
	noisemask = np.ones(len(thetTG))
	noisemask[(4,8,13),] = 0   # no red noise
	
	# Define a function that calculates the inverse of a singular vector decomposition
	def svdinv(traincv, y0):
		(m,n) = traincv.shape
		(U,S,V) = lin.svd(traincv, full_matrices=False)
		#s = np.diag(S)
		tol = np.max((m,n)) * np.abs(np.spacing(np.max(S)))
		r = np.sum(S > tol)
		#invtraincv = np.matmul(np.matmul(V[:,0:r], np.diag(S[0:r])**-1), U[:,0:r].T)
		invtraincv = np.matmul(np.matmul(V[:,0:r], S[0:r]**-1), U[:,0:r].T)
		alfa = np.matmul(invtraincv, y0)
		return(alfa, invtraincv, S, r)
	
	# Calculate the inverse singular vector decomposition
	print(traincv(thetTG))
	(alfa, invtraincv, s, r) = svdinv(traincv(thetTG), TGrsl)
	
	# Calculate the mean field
	fTG = testcv(thetTG*noisemask).T * alfa
	
	# Calculate V...whatever that is
	VTG = testcv2(thetTG*noisemask) - np.matmul(np.matmul(testcv(thetTG*noisemask), invtraincv), testcv(thetTG*noisemask).T)
	VTG = 0.5*(VTG+VTG.T)
	sdTG = np.sqrt(np.diag(VTG))
	
	return(TGtargcoords, fTG)