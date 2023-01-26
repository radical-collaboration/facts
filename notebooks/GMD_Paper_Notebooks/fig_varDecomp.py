import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
#
#
# ==================================================================================================================================
# Function to import nc data and preprocess for figure 6. 
#
def module_varFig_nc(df,WF,SCENARIO,EXPDATE):
    #  
    a = [];     components = []
    #
    for scenario in SCENARIO:
        for wf in WF: 
            #
            # Pick `Component/modules` for each workflow.
            COMPONENT = wf.split("-")[0]; MODULE = wf.split("-")[1]; 
            val = df.index[ (df['Component'] == COMPONENT) & (df['Module'] == MODULE) ].values[0]
            #       
            # Skip if data is absent.
            if df.at[val, 'DataFile'] in 'XXX':
                continue 
            #
            # Pick the data files & Import nc file to dataframe. 
            if SCENARIO[0][:3] == 'ssp':
                dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
                dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario) + df["DataFile"][val]
            elif SCENARIO[0][:3] == 'rcp':
                dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/ar5k14.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
                dataFILE    = 'ar5k14.{arg1}.'.format(arg1=scenario) + df["DataFile"][val]
            d_nc        = xr.open_dataset(dataFOLDER + dataFILE)
            #print( "DataFILE ==> "+str(df["DataFile"][val])+"\n"+ str(d_nc.dims) + "\n" + str(d_nc.data_vars)+ "\n" + str(d_nc.coords)+"\n"+"-------------"+"\n"+"\n" )
            #
            # Index for time.
            ST = 2020 ; EN = 2100
            YindST = np.where(d_nc["years"].values == ST)[0][0];   YindEN = np.where(d_nc["years"].values == EN)[0][0]
            # Save data into a new variable.
            components.append(COMPONENT)
            b = d_nc.sea_level_change[:,YindST:YindEN+1,0].values
            a.append(b[None,:] )
    #stack all at once
    sampsloccomponents = np.vstack(a);  sampsloccomponents = np.transpose(sampsloccomponents,(1,0,2))
    yrs=d_nc.years[YindST:YindEN+1].values
    #
    return sampsloccomponents, components, yrs; 
#
#
#
# ==================================================================================================================================
# Function to compute variance and fraction of variance. 
#
def varV_varF(sampsloccomponents,components):
    #
    varV = []; varF = [];  valIND = 0
    denom=np.var(np.sum(sampsloccomponents[:,:,:],axis = 1),axis=0)
    for co in components:
        # u = np.sum(sampsloccomponents[:,0:valIND+1,:],axis = 1)
        u = sampsloccomponents[:,valIND,:]
        VAR_V=(np.var(u,axis=0))/1e6
        VAR_F=(np.var(u,axis=0))/denom
        varV.append(VAR_V)
        varF.append(VAR_F)
        valIND += 1    
    return varV, varF