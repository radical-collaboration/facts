import numpy as np
import pandas as pd
import netCDF4 as nc
from pandas.plotting import table 
import xarray as xr
#
#
# ==================================================================================================================================
# Used to input module datasets for Figure 2 & 4
def module_Fig_nc(df,SCENARIO,EXPDATE,yrBREAK):
    #
    df_TEMP   = [] 
    #
    # Loop over .nc files.
    for val in df.index.values:
        for scenario in SCENARIO: 
            module      = df.at[val, 'Module']     
            sub_module  = df.at[val, 'subModule'] 
            component   = df.at[val, 'Component']
            # Skip if data is absent.
            if df.at[val, 'DataFile'] in 'XXX':
                continue  
            # Open the .nc data file. 
            dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
            dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario) + df["DataFile"][val]
            d_nc        = xr.open_dataset(dataFOLDER + dataFILE)
            #
            # Percentile calculation.
            percentList = [50, 5, 17, 83, 95]
            # Loop over years.
            for yy in d_nc["years"].values:
                if (yy > yrBREAK) or (yy ==2005):
                    continue
                else:
                    # Find year index to pick SLC value
                    yind = np.where(d_nc["years"].values == yy)[0][0]
                    # In meters.
                    GMSL = (d_nc["sea_level_change"][:,yind,:].values)/1000
                    # Find Percentile ranges.
                    pcntle = np.percentile(GMSL[:], percentList );    pcntle = np.around(pcntle,5)
                    #
                    df_TEMP.append( [component,module,sub_module, scenario, yy, ] + pcntle.tolist() )
    df_ptile = pd.DataFrame( df_TEMP, columns=['Component','Module','subModule', 'SSP', 'Year', ] + [ f'col_{x}' for xi, x in enumerate( percentList )] )
    return df_ptile
#
#
#
# ==================================================================================================================================
# Used to input module datasets for Figure 3 & 5
def wf_Fig_nc(df,SCENARIO,EXPDATE,yrBREAK):
    #
    df_TEMP   = [] 
    #
    # Loop over .nc files.
    for val in df.index.values:
        for scenario in SCENARIO: 
            component = df.at[val, 'Component']
            workflow  = df.at[val, 'Workflow']
            # Skip if data is absent.
            if df.at[val, 'DataFile'] in 'XXX':
                continue 
            # Skip for the following combinations. 
            if (scenario == 'ssp245') & (workflow == 'wf4'):
                continue
            if (scenario == 'ssp119') & ((workflow == 'wf3f') or (workflow == 'wf4')):
                continue
            if (scenario == 'ssp370') & ((workflow == 'wf3f') or (workflow == 'wf4')):
                continue
            #
            # Open the .nc data file. 
            dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
            dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario) + df["DataFile"][val]
            d_nc        = xr.open_dataset(dataFOLDER + dataFILE)
            #
            # Percentile calculation.
            percentList = [50, 5, 17, 83, 95]
            # Loop over years.
            for yy in d_nc["years"].values:
                if yy > yrBREAK:
                    continue
                else:
                    # Find year index to pick SLC value
                    Yind = np.where(d_nc["years"].values == yy)[0][0]
                    # In meters.
                    GMSL = (d_nc["sea_level_change"][:,Yind,:].values)/1000
                    # Find Percentile ranges.
                    pcntle = np.percentile(GMSL[:], percentList );    pcntle = np.around(pcntle,5)
                    #
                    df_TEMP.append( [component,workflow, scenario, yy, ] + pcntle.tolist() )
    df_ptile = pd.DataFrame( df_TEMP, columns=['Component','Workflow', 'SSP', 'Year', ] + [ f'col_{x}' for xi, x in enumerate( percentList )] )
    return df_ptile