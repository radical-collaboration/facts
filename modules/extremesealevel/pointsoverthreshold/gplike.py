#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" gplike.py
Python version of MATLAB's 'gplike':
    
    #GPLIKE Negative log-likelihood for the generalized Pareto distribution.
    #NLOGL = GPLIKE(PARAMS,DATA) returns the negative of the log-likelihood for
    #the two-parameter generalized Pareto (GP) distribution, evaluated at
    #parameters shape = K and scale = SIGMA, given DATA. GPLIKE does
    #not allow a threshold (location) parameter. NLOGL is a scalar.
    #
    #[NLOGL, ACOV] = GPLIKE(PARAMS,DATA) returns the inverse of Fisher's
    #information matrix, ACOV.  If the input parameter values in shape, scale are the
    #maximum likelihood estimates, the diagonal elements of ACOV are their
    #asymptotic variances.  ACOV is based on the observed Fisher's information,
    #not the expected information.
    #
    #When K = 0 and THETA = 0, the GP is equivalent to the exponential
    #distribution.  When K > 0 and THETA = SIGMA/K, the GP is equivalent to the
    #Pareto distribution.  The mean of the GP is not finite when K >= 1, and the
    #variance is not finite when K >= 1/2.  When K >= 0, the GP has positive
    #density for X>THETA, or, when K < 0, for 0 <= (X-THETA)/SIGMA <= -1/K.
    #
    #See also GPCDF, GPFIT, GPINV, GPPDF, GPRND, GPSTAT.
    
    #References:
    #      [1] Embrechts, P., C. Kl�ppelberg, and T. Mikosch (1997) Modelling
    #          Extremal Events for Insurance and Finance, Springer.
    #      [2] Kotz, S. and S. Nadarajah (2001) Extreme Value Distributions:
    #          Theory and Applications, World Scientific Publishing Company.
    
Parameters:
shape = shape coefficient of GPD
scale = scale coefficient of GPD
data  = data to which the GPD is fitted

Output: nlogl and acov, see above

Created on Tue Nov  5 12:24:26 2019
@author: Tim Hermans
"""
import numpy as np

def gplike(shape,scale,data):
    '''   

    '''
    k       = shape   # Tail index parameter
    sigma   = scale   # Scale parameter
    lnsigma = np.log(sigma);   # Scale parameter, logged
    
    n = len(data);
    z = data/sigma;
    
    if abs(k) > np.spacing(1):
        if k > 0 or max(z) < -1/k:
            u = 1 + k*z
            sumlnu = sum(np.log1p(k*z));
            nlogL = n*lnsigma + (1+1/k)*sumlnu;
            v = z/u;
            sumv = sum(v);
            sumvsq = sum(v**2);
            nH11 = 2*sumlnu/k**3 - 2*sumv/k**2 - (1+1/k)*sumvsq;
            nH12 = (-sumv + (k+1)*sumvsq)/sigma;
            nH22 = (-n + 2*(k+1)*sumv - k*(k+1)*sumvsq)/sigma**2;
            acov = [[nH22, -nH12],[-nH12, nH11]] / (nH11*nH22 - nH12*nH12)
        else:
            # The support of the GP when k<0 is 0 < y < abs(sigma/k)
            nlogL = np.Inf;
            acov = [[np.nan, np.nan], [np.nan, np.nan]];
            
    else: # limiting exponential dist'n as k->0
        # Handle limit explicitly to prevent (1/0) * log(1) == Inf*0 == NaN.
        nlogL = n*lnsigma + sum(z);
        sumz = sum(z);
        sumzsq = sum(z**2);
        sumzcb = sum(z**3);
        nH11 = (2/3)*sumzcb - sumzsq;
        nH12 = (-n + 2*sumz)/sigma**2;
        nH22 = (-sumz + sumzsq)/sigma;
        acov = [[nH22, -nH12], [-nH12, nH11]] / (nH11*nH22 - nH12*nH12)
    
    return nlogL,acov