import numpy as np
import pandas as pd
import netCDF4 as nc
from pandas.plotting import table 
import xarray as xr
#
#
# ==================================================================================================================================
# Used to input Workflow data.
#
def ncWF(DF,SCENARIO,EXPDATE):
    #
    modNO       = DF.index.values;   # Number of datasets Module outputs
    #
    a = [];     ccomp = []
    #        
    for scenario in SCENARIO:
        for val in modNO:
            WF = DF["Workflow"][val]
            # Skip if data is absent.
            if DF.at[val, 'DataFile'] in 'XXX':
                continue  
            #   
            # Pick the data files & Import nc file to dataframe. 
            dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
            dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario) + DF["DataFile"][val]
            d_nc        = xr.open_dataset(dataFOLDER + dataFILE)
            #
            # Index for time.
            ST = 2020 ; EN = 2100
            YindST = np.where(d_nc["years"].values == ST)[0][0];   YindEN = np.where(d_nc["years"].values == EN)[0][0]
            # Save data into a new variable.
            b = d_nc.sea_level_change[:,YindST:YindEN+1,0].values
            a.append(b[None,:] )
    #stack all at once
    SAMPSLOCRISE = np.vstack(a);  SAMPSLOCRISE = np.transpose(SAMPSLOCRISE,(1,2,0))
    yrs=d_nc.years[YindST:YindEN+1].values
    return SAMPSLOCRISE, yrs; 
    #
#
#
# ==================================================================================================================================
# GPDLogNExceedances
#
def GPDLogNExceedances(z,llambda,shape,scale,MHHW):
    #
    z0 = z
    z = np.maximum(z, 0)  #z[z<0] = 0
    #
    if shape<0:
        z=np.minimum(z,.99999*-scale/shape)
        logN = np.log(llambda*(1+(shape)*z/scale)**(-1/(shape)))
    elif shape==0:
        logN = np.log(llambda)-z/scale
    else:
        logN = np.log(llambda*(1+(shape)*z/scale)**(-1/(shape)))
    #    
    y= logN
    #
    #
    # for those points below threshold to MHHW, put on a Gumbel.
    # if np.max(MHHW.shape)>=1:
    # np.maximum(z0, MHHW(1), z0)
    z0=np.maximum(z0, MHHW) 
    sub = np.argwhere(z0<0)
    #if np.max(MHHW.shape)>=2:
    #    MHHWfreq=MHHW(2); 
    #else:
    MHHWfreq=365.25/2  
    # y(sub) = np.log(llambda)+(np.log(MHHWfreq)-np.log(llambda))*z0(sub)/MHHW(1)
    y[tuple(sub.T)] = np.log(llambda)+(np.log(MHHWfreq)-np.log(llambda))*z0[tuple(sub.T)]/MHHW
    #
    return y    
    #
#
#
# ==================================================================================================================================
# EFcurvWF
#
def EFcurvWF(PTILE,SAMPSLOCRISE,WFNO,yrs,testz,threshold,llambda,shape,scale):
    #
    workflo  = np.arange(WFNO)
    effcurve = np.empty((yrs.shape[0],testz.shape[0],WFNO,)); effcurve[:] = np.nan
    for wf in workflo: 
        sampslocrise = SAMPSLOCRISE[:,:,wf]
        #
        samps=[];
        samps = np.hstack((np.zeros((sampslocrise.shape[0],1)),sampslocrise)) / 1000
        #
        # samps=np.minimum(samps, np.quantile(samps, .999,axis=0)) #truncate samples viewed as physically implausible
        #
        #
        gg=[]; gg1=[];
        for tt in enumerate(yrs):
            # print(tt[0])
            gg  = testz[np.newaxis] - samps[:,tt[0]][np.newaxis].T
            if PTILE == 'mean': #Find mean
                gg1 = np.real(np.mean(np.exp(GPDLogNExceedances(gg-threshold,llambda,shape,scale,-threshold)),axis=0))
            else:               #Find Percentile
                gg1 = np.real(np.percentile(np.exp(GPDLogNExceedances(gg-threshold,llambda,shape,scale,-threshold)),PTILE,axis=0))
                #
            effcurve[tt[0],:,wf] = np.squeeze(gg1[np.newaxis].T)
            #
    return effcurve
    #
#
#
# ==================================================================================================================================
# 
#