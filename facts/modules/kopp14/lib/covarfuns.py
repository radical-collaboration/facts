import numpy as np
import numpy.matlib as npml
import scipy.special as sps

''' covarfuns.py

Multiple covariance functions to be used in Gaussian Process Regression

'''

def dYears(years1, years2):
	return(np.abs(np.subtract.outer(years1, years2)))

def dYears0(years1, years2):
	return(np.subtract.outer(years1, years2))

def angd(lat0, lon0, lat, lon):
	
	# Convert the input from degrees to radians
	(lat0, lon0, lat, lon) = np.deg2rad(np.array([lat0,lon0,lat,lon]))
	
	# Calculate the angle between the vectors
	temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + \
	(np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),\
	(np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))
	
	# Convert the results from radians to degrees and return
	return(np.rad2deg(temp))

def dDist(x1, x2, dmax=1000, fillval=1e6):
	temp = angd(np.transpose(npml.repmat(x1[:,0],x2.shape[0],1)),\
		np.transpose(npml.repmat(x1[:,1],x2.shape[0],1)),\
		npml.repmat(x2[:,0],x1.shape[0],1),\
		npml.repmat(x2[:,1],x1.shape[0],1))
	temp = temp + (fillval*(np.add.outer(x1[:,0], x2[:,0]) > dmax))
	return(temp)

def kMat1(dx, thetas):
	return(thetas[0]**2 * np.exp(-dx/thetas[1]))

def kMat3(dx, thetas):
	return(thetas[0]**2 * (1 + np.sqrt(3)*dx/thetas[1]) * np.exp(-np.sqrt(3)*dx/thetas[1]))

def kMat5(dx, thetas):
	return(thetas[0]**2 * (1 + (np.sqrt(5)*dx/thetas[1]) * (1+np.sqrt(5)*dx/thetas[1]/3)) * np.exp(-np.sqrt(5)*dx/thetas[1]))

def kSE(dx, thetas):
	return(thetas[0]**2 * np.exp(-(dx**2)/(2*thetas[1]**2)))

def kDELTA(dx, thetas):
	thetas = np.atleast_1d(thetas)
	return(thetas[0]**2 * (dx==0))

def kDP(years1, years2, thetas, GIAanchoryear=2005):
	thetas = np.atleast_1d(thetas)
	return(thetas[0]**2 * np.multiply.outer(years1-GIAanchoryear, years2-GIAanchoryear))

def kMatG(dx, thetas):
	eps = np.abs(np.spacing(1))
	fact1 = thetas[0]**2 * 2**(1-thetas[2])/sps.gamma(thetas[2])
	fact2 = (np.sqrt(2*thetas[2]) * (dx+eps)/thetas[1])**thetas[2]
	fact3 = sps.kn(thetas[2], np.sqrt(2*thetas[2])*(dx+eps)/thetas[1])
	return(fact1 * fact2 * fact3)

def kRQ(dx, thetas):
	return(thetas[0]**2 * (1+dx**2/(2*thetas[1]*thetas[2]))**(-thetas[2]))
	
def kDELTAG(ad, thetas):
	return(kDELTA(ad, thetas) * (ad<360))

def kGEOG(ad, thetas):
	return(kMat5(ad, thetas) * (ad<360))

#-----------------------------------------------------------------------------------------

def cvfuncL(t1,t2,dtlt2,thetas):
	return(kDP(t1,t2,thetas[0]) + kMatG(dtlt2,thetas[1:4]))

def cvfuncGLR(t1,t2,thetas,ad,dtlt2):
	term1 = cvfuncL(t1,t2,dtlt2,thetas[0:4])
	term2 = kDP(t1,t2,thetas[4]) * (kDELTAG(ad,np.sqrt(thetas[6])) + kGEOG(ad,np.append(np.sqrt(1-thetas[6]), thetas[5])))
	term3 = kMatG(dtlt2,thetas[7:10]) * (kDELTAG(ad,np.sqrt(thetas[10])) + kGEOG(ad,np.append(np.sqrt(1-thetas[10]), thetas[11])))
	term4 = kDELTAG(ad, np.sqrt(thetas[0]**2 + thetas[4]**2) * 50)
	return(term1 + term2 + term3 + term4)

def cvfuncS0(ad, thetas):
	return(thetas[0]**2 * (kDELTAG(ad,np.sqrt(thetas[2])) + kGEOG(ad, np.append(np.sqrt(1-thetas[2]), thetas[1]))))

def cvfuncS(ad, thetas):
	return(thetas[0]**2 * (kGEOG(ad, np.append(np.sqrt(thetas[2]), thetas[3])) + kGEOG(ad, np.append(np.sqrt(1-thetas[2]), thetas[1]))))