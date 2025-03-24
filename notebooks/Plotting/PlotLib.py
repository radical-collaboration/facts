import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

class PlotLib:

    def __init__(self,root_dir='/opt/facts', version_flag='Version Unspecified', figure_dim=[10,5]):
        self.version_flag = version_flag
        self.figure_dim = figure_dim

        # Sets up the necessary directories and paths
        self.root_dir = root_dir
        self.exp_dir = os.path.join(root_dir,'experiments')
        self.resource_dir = os.path.join(root_dir,'resources')
        self.out_dir = os.path.join(root_dir,'notebooks/Plotting/data')
        self.plot_dir = os.path.join(root_dir,'notebooks/Plotting/plots/')

        # Checks if the data directory is present in the plotting notebook directory and creates one if not
        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)

        # Checks if the plots directory is present in the plotting notebook directory and creates one if not
        if not os.path.isdir(self.plot_dir):
            os.mkdir(self.plot_dir)

        # Creates the datastructure to store the SSP quantile infromation
        ssp_data = {'module':[],'SSP-119':[],'SSP-126':[],'SSP-245':[],'SSP-370':[],'SSP-585':[]}
        temp_data = {'module':[],'1.5 C':[],'2.0 C':[],'3.0 C':[],'4.0 C':[],'5.0 C':[]}

        self.ssp_quantiles = pd.DataFrame(ssp_data)
        self.temp_quantiles = pd.DataFrame(temp_data)

        # Sets the datatype of the column to be strings
        for column in self.ssp_quantiles.columns:
            self.ssp_quantiles[column] = self.ssp_quantiles[column].astype('str')
        # Sets the datatype of the column to be strings
        for column in self.temp_quantiles.columns:
            self.temp_quantiles[column] = self.temp_quantiles[column].astype('str')

        self.import_modules_dict()



    def import_modules_dict(self):
         
        self.scenarios = ['119','126','245','370','585']
        self.module_dict = {
            'fair_gsat': ['temperature.fair.temperature_gsat.nc'],
            'emulandice_ais': ['emuAIS.emulandice.AIS_', 'EMULANDICE/AIS'],
            'emulandice_gris': ['emuGrIS.emulandice.GrIS_', 'EMULANDICE/GrIS'],
            'emulandice_glaciers': ['emuglaciers.emulandice.glaciers_',  'EMULANDICE/GLACIERS'],
            'emulandice2_ais': ['emuAIS.emulandice2.AIS_ALL_', 'EMULANDICE2/AIS'],
            'emulandice2_gris': ['emuGIS.emulandice2.GrIS_ALL_',  'EMULANDICE2/GrIS'],
            'emulandice2_glaciers': ['emuGLA.emulandice2.glaciers_ALL_',  'EMULANDICE2/GLACIERS'],
            'ipccar5_glaciers': ['ar5glaciers.ipccar5.glaciers_',  'IPCCAR5/GLACIERS'],
            'ipccar5_ais': ['ar5AIS.ipccar5.icesheets_AIS_',  'IPCCAR5/AIS'],
            'ipccar5_gris': ['ar5AIS.ipccar5.icesheets_GIS_',  'IPCCAR5/GrIS'],
            'larmip_ais': ['larmip.larmip.AIS_',  'LARMIP/AIS'],
            'fittedismip_gris': ['GrIS1f.FittedISMIP.GrIS_GIS_',  'FITTEDISMIP/GRIS'],
            'tlm_sterodynamics': ['ocean.tlm.sterodynamics_',  'TLM/STERODYNAMICS'],
            'bamber19_ais': ['bamber19.bamber19.icesheets_AIS_',  'BAMBER19/AIS'],
            'bamber19_gris': ['bamber19.bamber19.icesheets_GIS_',  'BAMBER19/GrIS'],
            'deconto21_ais': ['deconto21.deconto21.AIS_AIS_', 'DECONTO21/AIS'],
            'wf1e': ['total.workflow.wf1e.', 'WF1F'],
            'wf1f': ['total.workflow.wf1f.',  'WF1F'],
            'wf2e': ['total.workflow.wf2e.',  'WF2E'],
            'wf2f': ['total.workflow.wf2f.',  'WF2F'],
            'wf3e': ['total.workflow.wf3e.',  'WF3E'],
            'wf3f': ['total.workflow.wf3f.',  'WF3F'],
            'wf4': ['total.workflow.wf4.',  'WF4'],
            }  

    def get_module_data(self,filename, year=2100):
        # Function to obtain the GMSL projections for a specific module for a specified year as well as the quantiles
        # dat (string): the Path to the the data file wished to be open
        # year (int): The year to pull the data from

        data = (xr.open_dataset(filename).squeeze(drop=True).sea_level_change.sel(years=year, drop=True) / 1000).values
        module_data = []
        for i in range(len(data)):
            module_data.append(data[i])
        return module_data

   
    def get_gsat_data(self, filename, year=np.arange(2081, 2100)):
        # Function to obtain the GSAT surface air temperature from the FAIR Temperature Module
        # dat (string): the Path to the the data file wished to be open
        # year (int): The year to pull the data from

        data = (xr.open_dataset(filename).squeeze(drop=True).surface_temperature.sel(years=year, drop=True)).values
        gsat_data = []
        for i in range(len(data)):
            current_avg = np.average(data[i])
            gsat_data.append(current_avg)
        return gsat_data

    
    def plot_quantiles(self,temperatures, sea_levels, bin_start=1, bin_stop=7.0, bin_interval=0.5, cutoff=200,show=False):
        # Bins go from 1 C to 7.0 C in increments of 0.5 C
        # temperature (array): the GSAT data
        # sea_levels (array): The GMSL data
        # bin_start (float): Starting point for binning
        # bin_end (float): Ending point for binning
        # interval (float): Interval steps for binning
        # cutoff (int): The minimum number of samples for the bin to be plotted
        self.temp_quantiles.loc[self.mod_idx, 'module'] = self.module_dict[self.module][1]

        # Create bins based on the specified start, stop, and interval
        bin_centers = np.arange(bin_start, bin_stop, bin_interval)
        bins = [(center - bin_interval / 2) for center in bin_centers]
        bins.append(bins[-1] + bin_interval)

        # Digitize the data
        bin_idxs = np.digitize(temperatures, bins) - 1

        # Sets the color of the boxplot and the color of the median line within the box
        box_color = 'black'
        median_color = 'white'

        # Setting arbitrary high values for the ylimits so the script can correct them
        self.ylimits = [100,-100]
        table_centers = ['1.5', '2.0', '3.0', '4.0', '5.0' ]

        # Calculate and plot the quantiles for each bin
        for i, center in enumerate(bin_centers):
            # Get the sea levels for the current bin
            bin_data = [sea_levels[j] for j in range(len(sea_levels)) if bin_idxs[j] == i]

            if len(bin_data) >= cutoff:
                # Calculate the quantiles
                quantiles = self.get_quantiles(bin_data,show=False)
                current_quantiles = f'{np.round(quantiles[2],2)} ({np.round(quantiles[1],2)}-{np.round(quantiles[3],2)})'
                
                q5 = quantiles[0]
                q17 = quantiles[1]
                median = quantiles[2]
                q83 = quantiles[3]
                q95 = quantiles[4]

                if q5 <= self.ylimits[0]:
                    self.ylimits[0] = q5 - 0.02
                if q95 >= self.ylimits[1]:
                    self.ylimits[1] = q95 + 0.02

                # Plot the bar (median) and the whiskers (quantiles)
                plt.vlines(x=center, ymin=q5, ymax=q95, color='black')
                plt.boxplot(quantiles,
                            positions=[center],
                            sym="",
                            whis=0,
                            manage_ticks=False,
                            patch_artist=True,
                            boxprops=dict(facecolor=box_color, color=box_color),
                            capprops=dict(color=box_color),
                            whiskerprops=dict(color=box_color),
                            flierprops=dict(color=box_color, markeredgecolor=box_color),
                            medianprops=dict(color=median_color)
                            )
                
                if show:
                    print(f'{center}: {np.round(median,2)} ({np.round(q17,2)}-{np.round(q83,2)})')

                for table_center in table_centers:
                    if table_center == str(center):
                        self.temp_quantiles.loc[self.mod_idx, f'{center} C'] = current_quantiles

    
    # Gets the quantile information from inputted GMSL data
    def get_quantiles(self,data,label='',quantiles=[0.05, 0.17, 0.50, 0.83, 0.95],show=True):
        
        quantiles = np.quantile(data,quantiles)
        formatted_quants = [np.round(quantiles[2],2), np.round(quantiles[1],2), np.round(quantiles[3],2)]

        if show:
            print(f'{label.upper()}: {formatted_quants[0]} ({formatted_quants[1]}-{formatted_quants[2]})')
        
        return quantiles

    
    def plot_module(self, mod_idx, module, exp_name, localization='global'):
        
        # Sets the figure size as defined in the initialization function
        plt.figure(figsize=self.figure_dim)
        
        # Validate the module
        if module not in self.module_dict:
            raise ValueError(f"{module} is an invalid module name. Please choose from {list(self.module_dict.keys())}")

        self.module = module
        module_name = self.module_dict[self.module][1]
        
        self.mod_idx = mod_idx
        self.ssp_quantiles.loc[self.mod_idx, 'module'] = module_name
        
        # Plot settings
        xlim_range = [0.5, 5.5]

        plot_colors = {'119':'red', 
                       '126':'blue', 
                       '245':'green', 
                       '370':'orange', 
                       '585':'purple'}

        # Initialize lists to store the combined data
        combined_gsat = []
        combined_gmsl = []

        # Process each scenario
        for scenario in self.scenarios:
            this_scenario = f'{exp_name}{scenario}'
            
            # Checks to see if it's a workflow file as the paths do not inlclude the "sl" at the end 
            if "workflow" in self.module_dict[self.module][0]:
                self.module_path = f'{self.out_dir}/{this_scenario}.{self.module_dict[self.module][0]}{localization}.nc'
            else:
                self.module_path = f'{self.out_dir}/{this_scenario}.{self.module_dict[self.module][0]}{localization}sl.nc'

            # Get the GSAT and GMSL data
            gsat = self.get_gsat_data(f'{self.out_dir}/{this_scenario}.{self.module_dict["fair_gsat"][0]}')
            gmsl = self.get_module_data(filename=f'{self.module_path}')

            # Combine the data
            combined_gsat.extend(gsat)
            combined_gmsl.extend(gmsl)

            # Calculate the quantiles
            quants = self.get_quantiles(gmsl,label=f'{this_scenario}',show=False)
            current_quantiles = f'{np.round(quants[2],2)} ({np.round(quants[1],2)}-{np.round(quants[3],2)})'
            self.ssp_quantiles.loc[self.mod_idx, f'SSP-{scenario}'] = current_quantiles

            # Plot the data
            plt.scatter(x=gsat, 
                        y=gmsl, 
                        marker='o', 
                        s=40, 
                        color=plot_colors[scenario], 
                        alpha=0.1, 
                        edgecolors='none',
                        label=f'{scenario}: {current_quantiles}'
                        )

        # Overlay the binning data on the scatter plots
        self.plot_quantiles(
                            temperatures=combined_gsat, 
                            sea_levels=combined_gmsl,
                            show=False)

        # Set up the plot
        
        plt.xlim(*xlim_range)
        plt.ylim(self.ylimits[0], self.ylimits[1])
        plt.title(f'{module_name} {localization.upper()} [FACTS {self.version_flag}]\n NSAMPS PER SCENARIO = {2000}')
        plt.xlabel('2081-2100 Average GSAT [C$^\circ$]')
        plt.ylabel('2100 GMSL [m]')
        plt.legend(bbox_to_anchor=(1,1), loc='upper left')
        plt.tight_layout()

        # Save the plot
        plt.savefig(f'{self.plot_dir}/{module}_{localization}.png')
        plt.close()
