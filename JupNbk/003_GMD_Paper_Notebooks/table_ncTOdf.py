import numpy as np
import pandas as pd
import netCDF4 as nc
from pandas.plotting import table 
import xarray as xr
#
#
# ==================================================================================================================================
# Used to input module datasets for Table 3
def module_Table_nc(df,SCENARIO,SCENARIO_R,EXPDATE,YEAR):
    #
    df_TEMP   = [] 
    #
    # Loop over .nc files.
    for val in df.index.values:
        for scenario in enumerate(SCENARIO):
            #
            module     = df.at[val, 'Module']  
            sub_module = df.at[val, 'sub_Module']
            component  = df.at[val, 'Component'] 
            num        = df.at[val, 'Num']
            # Skip if data is absent.
            if df.at[val, 'DataFile'] in 'XXX':
                continue
            #
            # Open the .nc data file. 
            if df.at[val,'dataID'] == 'cplng':
                dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario[1],arg2=EXPDATE)
                dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario[1]) + df.at[val,'DataFile']
            elif df.at[val,'dataID'] == 'ar5k14':
                dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/ar5k14.{arg1}/output/'.format(arg1=SCENARIO_R[scenario[0]],arg2=EXPDATE)
                dataFILE    = 'ar5k14.{arg1}.'.format(arg1=SCENARIO_R[scenario[0]]) + df.at[val,'DataFile']
                #
            d_nc = xr.open_dataset(dataFOLDER + dataFILE)
            #
            # Percentile calculation.
            percentList = [50, 17, 83]      
            # Find year index to pick SLC value
            if sub_module == 'thermalexpansion':
                Yind = np.where(d_nc["year"].values == YEAR)[0][0]
            else:
                Yind = np.where(d_nc["years"].values == YEAR)[0][0]
            # 
            if sub_module == 'temperature':
                GMSL = (d_nc["surface_temperature"][:,Yind,0].values) 
            else:
                GMSL = (d_nc["sea_level_change"][:,Yind,0].values)/1000       
            # Find Percentiles.
            pcntle = np.percentile(GMSL[:], percentList );    pcntle = np.around(pcntle,2)
            #
            MODULE_PAPER   = module+'/'+sub_module
            SCENARIO_PAPER = scenario[1]+'/'+SCENARIO_R[scenario[0]]
            df_TEMP.append( [num,component,MODULE_PAPER, SCENARIO_PAPER] + pcntle.tolist() )
            #
    df_ptile = pd.DataFrame( df_TEMP, columns=['Num','Component','Module','SSP/RCP', ] + [ f'col_{x}' for xi, x in enumerate( percentList )] )
    df_ptile[''] = df_ptile.apply(lambda x: f'{x.col_50:2.2f} ({x.col_17:2.2f} - {x.col_83:2.2f})', axis=1 )
    df_ptile1    = pd.DataFrame( df_ptile.set_index( ['Num','Component','Module','SSP/RCP'] )[''] ).unstack().swaplevel( 0,1, axis=1 )
    return df_ptile1
#
#
#
# ============================================================================================================================================================
# Used to input workflow datasets for Table 3
def wf_Table_nc(df,SCENARIO,EXPDATE,YEAR):
    #
    df_TEMP   = [] 
    #
    # Loop over .nc files.
    for val in df.index.values: 
        for scenario in SCENARIO: 
            workflow  = df.at[val, 'Workflow'];  
            component = df.at[val, 'Component'];  
            # Skip if data is absent.
            if df.at[val, 'DataFile'] in 'XXX':
                continue 
            if (scenario == 'ssp119') & ((workflow == 'wf3e') or (workflow == 'wf3f') or (workflow == 'wf4')):
                continue
            if (scenario == 'ssp245') & ((workflow == 'wf4')):
                continue
            if (scenario == 'ssp370') & ((workflow == 'wf3e') or (workflow == 'wf3f') or (workflow == 'wf4')):
                continue
            # Open the .nc data file.
            dataFOLDER  = '/projects/kopp/facts-experiments/{arg2}/coupling.{arg1}/output/'.format(arg1=scenario,arg2=EXPDATE)
            dataFILE    = 'coupling.{arg1}.'.format(arg1=scenario) + df.at[val, 'DataFile']
            d_nc        = xr.open_dataset(dataFOLDER + dataFILE)
            #
            # Percentile calculation.
            percentList = [50, 17, 83]        
            # Find year index to pick SLC value
            Yind    = np.where(d_nc["years"].values == YEAR)[0][0]
            GMSL    = (d_nc["sea_level_change"][:,Yind,0].values)/1000  # Convert to meters.
            # Find Percentiles.
            pcntle  = np.percentile(GMSL[:], percentList );    pcntle = np.around(pcntle,2)
            #
            df_TEMP.append( [component,workflow, scenario] + pcntle.tolist() )
    df_ptile = pd.DataFrame( df_TEMP, columns=['Component','Workflow','SSP', ] + [ f'col_{x}' for xi, x in enumerate( percentList )] )
    #
    df_ptile[''] = df_ptile.apply(lambda x: f'{x.col_50:2.2f} ({x.col_17:2.2f} - {x.col_83:2.2f})', axis=1 )
    df_ptile1    = pd.DataFrame( df_ptile.set_index( ['Component', 'Workflow', 'SSP'] )[''] ).unstack().swaplevel( 0,1, axis=1 )
    return df_ptile1