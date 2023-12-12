import pickle
import sys
import numpy as np
from scipy.stats import norm

from cholcov import cholcov

''' SampleISDists.py

This function produces samples from the fitted ice sheet distributions produced in the
fitting/calibration stage.

Parameters: 
nsamps = Number of samples to produce
sigmas = 1D array of sigmas for each ice sheet component
mus = 1D array of mus for each ice sheet component
offsets = 1D array of offsets for each ice sheet component
islastdecade = Ice sheet accelerations(?) for the last decade
corris = Correlation structure across ice sheets
seed = Seed for the random number generator

Return: An array of dimensions [nsamps, length(sigmas)] of samples of ice sheet melts
accelerations.

Note: The lengths of 'sigmas', 'mus', 'offsets', and 'islastdecade' must be the same.
Also, 'corris' must be a square matrix of the same size as the length of the previous 
parameters. 

'''

def SampleISDists(nsamps, sigmas, mus, offsets, islastdecade, corris, seed=1234):
	
	# Evenly sample an inverse normal distribution
	rng = np.random.default_rng(seed)
	x = np.linspace(0,1,nsamps+2)[1:(nsamps+1)]
	norm_inv = norm.ppf(x)
	
	# Find the Cholesky Covariance
	covis = np.dot(np.dot(np.diag(sigmas), corris), np.diag(sigmas))
	T = cholcov(covis)
	
	# Build a matrix of permutated norm_inv values
	norm_inv_perm = np.full((nsamps, sigmas.shape[0]), np.nan)
	for i in range(sigmas.shape[0]):
		norm_inv_perm[:,i] = rng.permutation(norm_inv)
	
	# Create the correlated samples for the ice sheet accelerations
	sampeps = np.dot(norm_inv_perm, T)
	sampisrates = offsets + np.exp((sampeps + mus))
	sampisaccel = (sampisrates - islastdecade) / (2100 - 2011)
	
	# Return the sampled accelerations
	return(sampisaccel)