import numpy as np

##########################################################################################
#
# sample_from_quantiles
#
# Produces a 1-d array of samples from a distribution represented by quantile values using
# linear interpolation between quantile values.
#
# Parameters:
#   qvals  - The quantile values that represent the distribution
#   q      - The quantiles at which 'qvals' are determined
#   nsamps - Number of samples to draw from this distribution
#
# Note: This sampling approach works better the more quantiles you provide, especially
# out in the tails of the distribution.
#
##########################################################################################
def sample_from_quantiles(qvals, q, nsamps):
	
	pool = np.interp(np.linspace(0,1,nsamps), q, qvals)
	rng.shuffle(pool)
	return(pool)
