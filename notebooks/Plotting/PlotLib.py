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
        self.temp_qauntiles = pd.DataFrame(temp_data)

        # Sets the datatype of the column to be strings
        for column in self.ssp_quantiles.columns:
            self.ssp_quantiles[column] = self.ssp_quantiles[column].astype('str')
        # Sets the datatype of the column to be strings
        for column in self.temp_qauntiles.columns:
            self.temp_qauntiles[column] = self.temp_qauntiles[column].astype('str')

        self.import_modules_dict()



    def import_modules_dict(self):
         
        self.scenarios = ['119','126','245','370','585']
        
        # Note for mod_names calls
        # mod_names['mod_name'][0] -> file tail for the module provided
        # mod_names['mod_name'][1] -> plot lower limit
        # mod_names['mod_name'][2] -> plot upper limit
        # mod_names['mod_name'][3] -> Title  
        self.module_dict = {
            'fair_gsat': ['temperature.fair.temperature_gsat.nc'],
            'emulandice_ais': ['emuAIS.emulandice.AIS_globalsl.nc', -0.05, .25, 'EMULANDICE/AIS'],
            'emulandice_gris': ['emuGrIS.emulandice.GrIS_globalsl.nc', -0.05, .225, 'EMULANDICE/GrIS'],
            'emulandice_glaciers': ['emuglaciers.emulandice.glaciers_globalsl.nc', .03, .23, 'EMULANDICE/GLACIERS'],
            'emulandice2_ais': ['emuAIS.emulandice2.AIS_ALL_globalsl.nc', -4.5, 3.5, 'EMULANDICE2/AIS'],
            'emulandice2_gris': ['emuGIS.emulandice2.GrIS_ALL_globalsl.nc', -2, 1.5, 'EMULANDICE2/GrIS'],
            'emulandice2_glaciers': ['emuGLA.emulandice2.glaciers_ALL_globalsl.nc', -.75, .65, 'EMULANDICE2/GLACIERS'],
            'ipccar5_glaciers': ['ar5glaciers.ipccar5.glaciers_globalsl.nc', .02, .275, 'IPCCAR5/GLACIERS'],
            'ipccar5_ais': ['ar5AIS.ipccar5.icesheets_AIS_globalsl.nc', -.1, .2, 'IPCCAR5/AIS'],
            'ipccar5_gris': ['ar5AIS.ipccar5.icesheets_GIS_globalsl.nc', .01, .35, 'IPCCAR5/GrIS'],
            'larmip_ais': ['larmip.larmip.AIS_globalsl.nc', -.05, .65, 'LARMIP/AIS'],
            'fittedismip_gris': ['GrIS1f.FittedISMIP.GrIS_GIS_globalsl.nc', 0, .25, 'FITTEDISMIP/GRIS'],
            'tlm_sterodynamics': ['ocean.tlm.sterodynamics_globalsl.nc', .07, .45, 'TLM/STERODYNAMICS'],
            'bamber19_ais': ['bamber19.bamber19.icesheets_AIS_globalsl.nc', -.20, 1.50, 'BAMBER19/AIS'],
            'bamber19_gris': ['bamber19.bamber19.icesheets_GIS_globalsl.nc', 0, 1.03, 'BAMBER19/GrIS'],
            'deconto21_ais': ['deconto21.deconto21.AIS_AIS_globalsl.nc', .05, .65, 'DECONTO21/AIS'],
            'wf1e_global': ['total.workflow.wf1e.global.nc', 0.15, 0.95, 'WF1F/GLOBAL'],
            'wf1f_global': ['total.workflow.wf1f.global.nc', 0.15, 0.95, 'WF1F/GLOBAL'],
            'wf2e_global': ['total.workflow.wf2e.global.nc', 0.15, 0.95, 'WF2E/GLOBAL'],
            'wf2f_global': ['total.workflow.wf2f.global.nc', 0.15, 0.95, 'WF2F/GLOBAL'],
            'wf3e_global': ['total.workflow.wf3e.global.nc', 0.15, 0.95, 'WF3E/GLOBAL'],
            'wf3f_global': ['total.workflow.wf3f.global.nc', 0.15, 0.95, 'WF3F/GLOBAL'],
            'wf4_global': ['total.workflow.wf4.global.nc', 0.15, 0.95, 'WF4/GLOBAL'],
            'wf1e_local': ['total.workflow.wf1e.local.nc', 0.15, 0.95, 'WF1F/LOCAL'],
            'wf1f_local': ['total.workflow.wf1f.local.nc', 0.15, 0.95, 'WF1F/LOCAL'],
            'wf2e_local': ['total.workflow.wf2e.local.nc', 0.15, 0.95, 'WF2E/LOCAL'],
            'wf2f_local': ['total.workflow.wf2f.local.nc', 0.15, 0.95, 'WF2F/LOCAL'],
            'wf3e_local': ['total.workflow.wf3e.local.nc', 0.15, 0.95, 'WF3E/LOCAL'],
            'wf3f_local': ['total.workflow.wf3f.local.nc', 0.15, 0.95, 'WF3F/LOCAL'],
            'wf4_local': ['total.workflow.wf4.local.nc', 0.15, 0.95, 'WF4/LOCAL'],
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

    
    def plot_quantiles(self,module,temperatures, sea_levels, bin_start=1, bin_stop=7.0, bin_interval=0.5, cutoff=200,show=False):
        # Bins go from 1 C to 7.0 C in increments of 0.5 C
        # temperature (array): the GSAT data
        # sea_levels (array): The GMSL data
        # bin_start (float): Starting point for binning
        # bin_end (float): Ending point for binning
        # interval (float): Interval steps for binning
        # cutoff (int): The minimum number of samples for the bin to be plotted

        #
        self.temp_qauntiles.loc[self.mod_idx, 'module'] = self.module_dict[self.module][3]

        # Create bins based on the specified start, stop, and interval
        bin_centers = np.arange(bin_start, bin_stop, bin_interval)
        bins = [(center - bin_interval / 2) for center in bin_centers]
        bins.append(bins[-1] + bin_interval)

        # Digitize the data
        bin_idxs = np.digitize(temperatures, bins) - 1

        # Sets the color of the boxplot and the color of the median line within the box
        box_color = 'black'
        median_color = 'white'

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
                        self.temp_qauntiles.loc[self.mod_idx, f'{center} C'] = current_quantiles

    
    # Gets the quantile information from inputted GMSL data
    def get_quantiles(self,data,label='',quantiles=[0.05, 0.17, 0.50, 0.83, 0.95],show=True):
        
        quantiles = np.quantile(data,quantiles)
        formatted_qaunts = [np.round(quantiles[2],2), np.round(quantiles[1],2), np.round(quantiles[3],2)]

        if show:
            print(f'{label.upper()}: {formatted_qaunts[0]} ({formatted_qaunts[1]}-{formatted_qaunts[2]})')
        
        return quantiles

    
    def plot_module(self, mod_idx, module, exp_name, use_ssp_tag=True):
        
        # Sets the figure size as defined in the initialization function
        plt.figure(figsize=self.figure_dim)
        
        # Validate the module
        if module not in self.module_dict:
            raise ValueError(f"{module} is an invalid module name. Please choose from {list(self.module_dict.keys())}")

        self.module = module
        module_name = self.module_dict[self.module]
        self.mod_idx = mod_idx
        self.ssp_quantiles.loc[self.mod_idx, 'module'] = module_name[3]
        
        # Plot settings
        plot_title = module_name[3]
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

            # Get the GSAT and GMSL data
            gsat = self.get_gsat_data(f'{self.out_dir}/{this_scenario}.{self.module_dict["fair_gsat"][0]}')
            gmsl = self.get_module_data(filename=f'{self.out_dir}/{this_scenario}.{module_name[0]}')

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
        self.plot_quantiles(module=module_name,
                            temperatures=combined_gsat, 
                            sea_levels=combined_gmsl,
                            show=False)

        # Set up the plot
        
        plt.xlim(*xlim_range)
        plt.ylim(self.ylimits[0], self.ylimits[1])
        plt.title(f'{plot_title} {self.version_flag}\n NSAMPS PER SCENARIO = {2000}')
        plt.xlabel('2081-2100 Average GSAT [C$^\circ$]')
        plt.ylabel('2100 GMSL [m]')
        plt.legend(bbox_to_anchor=(1,1), loc='upper left')
        plt.tight_layout()

        # Save the plot
        plt.savefig(f'{self.plot_dir}/{module}.png')
        plt.close()

