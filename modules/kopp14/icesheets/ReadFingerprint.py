import numpy as np
import sys
import pyssht as ssht

''' ReadFingerprint.py

Provides a function that reads in a fingerprint data file and converts the sine/cosine
coefficients to a signal on a lat/lon grid

Parameters:
fname = File name that contains the sine/cosine coefficients for the fingerprint

Return:
f = The fingerprint coefficient along a lat/lon grid [nlat x nlon]
lat = Vector of latitudes
lon = Vector of longitudes

Note: This function requires an external library "pyssht". Source code and instructions
for compilation and installation can be found at the url below:

http://astro-informatics.github.io/ssht/

Subsequently, the python library must be contained in your $PYTHONPATH.

'''

def ReadFingerprint(fname):

	# Calculates the number of real spherical harmonic orders that belong to an
	# expansion from degree l=0 to L.
	def addmup(x):
		x = np.array(x)
		res = 0.5*(x+1)**2+0.5*x+0.5
		return(res.astype(int))
	
	# Read in the file header 
	with open(fname) as f:
		header = f.readline()
		(L, age) = header.split()
		L = int(L)

	# Read in the rest of the file
	test = np.loadtxt(fname, skiprows=1)
	my_C = test[:,[0,2]].flatten()
	my_S = test[:,[1,3]].flatten()

	# Calculate the indices that are expected to have '0' as coefficients
	mods = L - np.arange(0,L+1) + 1
	modu = np.cumsum(mods + np.mod(mods, 2))
	zero_inds = modu[np.mod(mods,2) == 1] - 1

	# Remove the zero coefficients
	my_C = np.delete(my_C, zero_inds)
	my_S = np.delete(my_S, zero_inds)

	# Indices to line up coefficients in hlm
	modm = np.append(np.cumsum(np.append(0, mods[np.arange(0,mods.size-1)])), int(addmup(L)))

	#### Produces "dems' and 'dels' from addmon
	dems = np.empty(0, np.int8)
	dels = np.empty(0, np.int8)

	for i in np.arange(0,L+1):
		dems = np.append(dems, np.arange(0,i+1))
		dels = np.append(dels, np.repeat(i, i+1))
	
	# Initialize variable to hold coefficients
	hlm = np.zeros((dels.size,2))

	# Populate hlm with the coefficients
	for m in np.arange(0,L+1):
		hlm[addmup(np.arange(m-1, L))+m,0] = my_C[modm[m]:modm[m+1]]
		hlm[addmup(np.arange(m-1, L))+m,1] = my_S[modm[m]:modm[m+1]]

	# Realign conventions
	CSC = (-1)**dems
	hlm[:,0] = hlm[:,0]*CSC
	hlm[:,1] = hlm[:,1]*CSC

	# Convert hlm into a form that can be used by the ssht library
	counter = 0
	flm = np.zeros((L+1)**2, dtype=np.complex_)
	for ii in np.arange(0,L+1):
		sub = np.where(dels == ii)[0]
		if sub.size > 1:
			sub2 = np.append(sub[-1:0:-1], sub)
		else:
			sub2 = sub
		flm[counter:(counter+sub2.size)] = hlm[sub2,0] + 1j*hlm[sub2,1]
		counter = counter + sub2.size

	fels = np.floor(np.sqrt(np.arange(0,flm.size)))
	fems = np.arange(0,flm.size) - fels*fels - fels
	sub = np.where(fems < 0)[0]
	flm[sub] = (-1)**fems[sub] * np.conj(flm[sub])

	# Invert the harmonics
	f = ssht.inverse(flm, L+1, Method="MWSS", Reality=True)
	f = f*np.sqrt(4*np.pi)
	(thetas, phis) = ssht.sample_positions(L+1, "MWSS")
	lat = 90-np.rad2deg(thetas)
	lon = np.mod(np.rad2deg(phis)+180, 360)
	sort_inds = np.argsort(lon)
	lon = lon[sort_inds]
	f = f[:,sort_inds]

	return(f, lat, lon)