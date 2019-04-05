import numpy as np
import scipy.linalg as sp

def cholcov(x):
	
	# Determine the rank of the matrix
	rank = np.linalg.matrix_rank(x)
	
	# If the rank of the matrix is not equal to the number of rows/columns of the matrix
	if rank != x.shape[0]:
		
		# Calculate the eigen values/vectors
		vals, vecs = sp.eig(x)
		
		# Keep vectors whose values are greater than zero (within tolerance)
		red_vals = np.sqrt(vals[vals>1e-8])
		red_vecs = vecs[:,vals>1e-8]
		
		# Calculate T
		T = np.dot(np.diag(red_vals), np.matrix.transpose(red_vecs))
	
	else:
		
		# Calculate the cholesky factors
		T = sp.cholesky(x)
	
	return(T)