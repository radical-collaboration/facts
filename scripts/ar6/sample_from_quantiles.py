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

	pool = np.interp(np.linspace(0,1,nsamps+2)[1:(nsamps+1)], q, qvals)
	np.random.shuffle(pool)
	return(pool)


if __name__ == "__main__":

	# Generate random data
	rng = np.random.default_rng()
	x = rng.normal(10, 2, (120, 1000))	# 120 "times", 1000 "samples"
	x += np.linspace(0,59, 120)[:,None]

	qs = np.linspace(0,1,102)[1:101]
	qvals = np.quantile(x, qs, axis=1)

	temp = []
	np.random.seed(1234)
	for i in np.arange(120):
		#np.random.seed(1234)	# Uncomment to get correlated over "time" dimension
		temp.append(sample_from_quantiles(qvals[:,i], qs, 100))

	temp = np.array(temp)
	print(np.diff(temp[:,0]))
