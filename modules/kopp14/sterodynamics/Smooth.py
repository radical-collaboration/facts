import numpy as np

''' Smooth.py

Replicates MatLab's "smooth" function.

Parameters: 
x = Vector of data to smooth
w = Number of years contained in the smoothing window (Must be odd)

Return: 
y = Smoothed vector

'''

def Smooth(x, w=5):
	out0 = np.convolve(x, np.ones(w,dtype='double'), 'valid')/w
	r = np.arange(1,w-1,2, dtype="double")
	start = np.cumsum(x[:w-1])[::2]/r
	stop = (np.cumsum(x[:-w:-1])[::2]/r)[::-1]
	y = np.concatenate((start, out0, stop))
	return(y)