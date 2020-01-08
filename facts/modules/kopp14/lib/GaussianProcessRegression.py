import numpy as np
import numpy.linalg as lin

''' GaussianProcessRegression.py

This script applies Gaussian Process Regression.  This is adapted from the matlab
code of the same name in Kopp 2104 workflow, which was "modeled" on Rasmussen &
Williams 2006

Parameters: 
x0 = x values at training points
y0 = y values at training points (demeaned)
x = x values at which to approximate
traincv = covariance matrix among covariance points
cvfunc = EITHER the handle to a function that calculates covariance between x0 and x
		 or the covariance matrix between x0 and x
testcv2 = If cvfunc is passed as a matrix, then the covariance matrix among x
invcv = "structure with invcv.alfa, invcv.invtraincv, invcv.s, invcv.r"...whatever that
		means.

Return: 
f = 
V = 
logp = 
alfa = 
errorflags = 
invtraincv = 
invcv = 

Note: 
No description was provided for any of the return variables. 

'''

def GaussianProcessRegression(x0, y0, x, traincv, cvfunc, testcv2, invcv):
	
	
	
	return()