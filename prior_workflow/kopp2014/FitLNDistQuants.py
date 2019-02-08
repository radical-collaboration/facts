import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

'''
FitLNDistQuants

Fit parameters of a log-normal distribution using the provided quantile values.

Parameters:
target_mean_val = The initial guess of the distribution mean
target_quant_vals = The target values for the fitted the 'quant_levs'.
minval = A shift in the mean used in the initial guess
quant_levs = A list of the desired quantile levels for the fitting process

Return:
thet = The parameters of the fitted log-normal distribution

Note: This function is transcoded from Robert Kopp's FitLNDistributionQuantiles.m. This
also assumes that target_quant_vals[1] is the 50th percentile value so as to estimate an
initial guess for mu.

'''
  
def FitLNDistQuants(target_mean_val, target_quant_vals, minval, quant_levs):
	
	# Returns the quantile values for the proposed set of distribution parameters
	def targetquants(x0, mu, sigma):
		return np.array(x0 + np.exp(mu + norm.ppf(quant_levs) * sigma))
	
	# Come up with an initial guess for the mean of the distribution
	muguess = np.log(target_quant_vals[1] - minval)
	
	# Come up with an initial guess for sigma
	if(target_mean_val <= target_quant_vals[1]):
		sigmaguess = 1.0
	else:
		sigmaguess = np.sqrt(np.log((target_mean_val - minval) / np.exp(muguess))*2)
	
	# Minimize the sum-of-squared residuals to estimate the distribution parameters
	def mintarget(thet):
		return np.sum((targetquants(thet[0], thet[1], thet[2]) - target_quant_vals)**2)
	
	# Create a vector of the initial guess
	init_guess = [minval, muguess, sigmaguess]
	
	# Define the sampling bounds in the optimization process
	my_bnds = ((-100E6,100E6),(-100E6, 100E6), (1E-9, 100))

	# Run the global optimization to estimate the parameters
	optim_results = minimize(mintarget, method="L-BFGS-B", x0=init_guess, bounds=my_bnds)
	
	# Return the parameter values
	return optim_results.x


