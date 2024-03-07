import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4


class FactsUtils:
    def __init__(self) -> None:
        ssp119_ref = {'bamber19_ais': [-0.01, 0.10, 0.27],
                      'bamber19_gis': [0.07, 0.13, 0.32],
                      'deconto21_ais': [0.07, 0.09, 0.12],
                      }
        
        ssp126_ref = {'bamber19_ais': [-0.01, 0.11, 0.31],
                      'bamber19_gis': [0.07, 0.14, 0.35],
                      'deconto21_ais': [0.07, 0.09, 0.12],
                      }
        
        ssp245_ref = {'bamber19_ais': [-0.00, 0.14, 0.43],
                      'bamber19_gis': [0.08, 0.16, 0.44],
                      'deconto21_ais': [0.07, 0.10, 0.29],
                      }
        
        ssp370_ref = {'bamber19_ais': [0.01, 0.17, 0.5],
                      'bamber19_gis': [0.09, 0.19, 0.53],
                      'deconto21_ais': [0.09, 0.22, 0.47],
                      }
        
        ssp585_ref = {'bamber19_ais': [0.02, 0.19, 0.55],
                      'bamber19_gis': [0.09, 0.21, 0.57],
                      'deconto21_ais': [0.12, 0.30, 0.50],
                      }

        self.reference_dicts = {'ssp119': ssp119_ref,
                                'ssp126': ssp126_ref,
                                'ssp245': ssp245_ref,
                                'ssp370': ssp370_ref,
                                'ssp585': ssp585_ref
                                }
        

    # This function should be placed in the projection stage of the module.
    # Note this automatically converts the runs projection data to m
    def verify_module_output(self, scenario, module_set, region, filename, year=2100):

        quantiles = [.17, .50, .83]

        # Opens the projection.nc file for the run data
        run_data = xr.open_dataset(filename)

        # loads the reference dictionary for this SSP
        reference_dict = self.reference_dicts[f'{scenario}'][f'{module_set}_{region}']

        # Create the pandas dataframe to output quantile information
        quantile_data = {'RUN (m)': [], 'REF (m)': [], 'DIFF (%)': []}
        quantile_data = pd.DataFrame(quantile_data)

        # Loops through the quantiles and adds them to the data frame
        for i in range(len(quantiles)):
            ref_quant = reference_dict[i]
            run_quant = self.get_quantile(data=run_data, year=year, quantile=quantiles[i])
            diff = np.round(((run_quant - ref_quant) / run_quant) * 100,2)
            quantile_data.loc[i, 'RUN (m)'] = run_quant
            quantile_data.loc[i, 'REF (m)'] = ref_quant
            quantile_data.loc[i, 'DIFF (%)'] = diff

        # Renames the indicies of the dataframe
        quantile_data.index = ['Q17', 'Q50', 'Q83']

        # SAves the check to the tasks directory int he sandbox
        quantile_data.to_csv(f'{module_set}_{region}.txt', sep=' ')
        return(quantile_data)

    def get_quantile(self, data, year, quantile, conversion=1000, round_quant=True):
        quantile = float(
            data.sea_level_change.sel(years=year, drop=True).quantile([quantile], dim='samples') / conversion)
        if round_quant:
            quantile = np.round(quantile, 2)
        return quantile
        