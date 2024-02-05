import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

class PlotLib:

    def __init__(self,root_dir='/opt/facts', version_flag='Version Unspecified'):
        self.version_flag = version_flag

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
        
       

    # Function to obtain the GMSL projections for a specific module for a specified year as well as the quantiles
    # dat (string): the Path to the the data file wished to be open
    # year (int): The year to pull the data from
    def get_module_data(self,filename, year=2100):
        data = (xr.open_dataset(filename).squeeze(drop=True).sea_level_change.sel(years=year, drop=True) / 1000).values
        module_data = []
        for i in range(len(data)):
            module_data.append(data[i])
        return module_data

    # Function to obtain the GSAT surface air temperature from the FAIR Temperature Module
    # dat (string): the Path to the the data file wished to be open
    # year (int): The year to pull the data from
    def get_gsat_data(self, filename, year=np.arange(2081, 2100)):
        data = (xr.open_dataset(filename).squeeze(drop=True).surface_temperature.sel(years=year, drop=True)).values
        gsat_data = []
        for i in range(len(data)):
            current_avg = np.average(data[i])
            gsat_data.append(current_avg)
        return gsat_data


    # Bins go from 0.0 C to 8.0 C in increments of 0.5 C
    # data1 (array): the GSAT data
    # data2 (array): The GMSL data
    # bin_start (float): Starting point for binning
    # bin_end (float): Ending point for binning
    # interval (float): Interval steps for binning
    # cutoff (int):
    # term_out (bool): prints outputs to terminal
    # plot_ax (bool): Plots vertical ax line at center of bin for diagnostics
    def bin_data(self, data1, data2, bin_start=0.25, bin_stop=5, interval=0.5, cutoff=5, term_out=False, plot_ax=True):

        # bins the data and creates a pandas dataframe
        bins = np.arange(start=bin_start, stop=bin_stop, step=interval)
        bin_idxs = np.arange(start=0, stop=len(bins) - 1)

        # Creates a Pandas Dataframe and defines the values to go in the individual bins outlined above
        dataframe = pd.DataFrame({"gsat": data1, "gmsl": data2})
        dataframe['binID'] = pd.cut(dataframe['gsat'], bins, include_lowest=True, labels=bin_idxs)

        for bin_idx in range(len(bin_idxs)):
            counter = 0
            for gmsl_idx in range(len(dataframe['gsat'])):
                if dataframe['binID'][gmsl_idx] == bin_idx:
                    counter += 1

        print('QUANTILES BY TEMPERATURE BIN')
        for bin_idx in range(len(bin_idxs)):
            current_bin = []
            current_pos = bins[bin_idx] + (interval / 2)

            # Appends the current bin with all the values from the dataframe for _this_ bin
            for gmsl_idx in range(len(dataframe['gmsl'])):
                if dataframe['binID'][gmsl_idx] == bin_idx:
                    current_bin.append(dataframe['gmsl'][gmsl_idx])

            if len(current_bin) >= cutoff:
                # Gets the proper quantiles from np.quantile
                bin_quants = np.quantile(current_bin, [.05, .17, .5, .83, .95])
                print(f"{bin_idx} {current_pos}C: {np.round(bin_quants[2],2)} ({np.round(bin_quants[1],2)}-{np.round(bin_quants[3],2)})")
                box_color = 'black'
                median_color = 'white'

                for i in range(len(bin_quants)):
                    plt.vlines(x=current_pos, ymin=bin_quants[0], ymax=bin_quants[4], color='black')

                # Plots a box whisker plot of the quantiles defined above
                plt.boxplot(bin_quants,
                            positions=[current_pos],
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

            if plot_ax:
                for i in range(len(bin_idxs)):
                    plt.axvline(bins[i] + (interval / 2),
                                linewidth=0.5,
                                linestyle="--",
                                color=(0, 0, 0, 0.01))

        return dataframe

    # Gets the quantile information from inputted GMSL data
    def get_quants(self, data, ssp_label='',save_quants=True):
        quants = np.quantile(data, [.05, .17, .5, .83, .95])
        quants_out = f'{np.round(quants[2],2)} ({np.round(quants[1],2)}-{np.round(quants[3],2)})'
        print(f'{ssp_label.upper()}: {quants_out}')

        return quants_out

    
    def plot_module(self, module, exp_name, use_ssp_tag=True):
        # Checks to see if module exists in the module dictionary
        if module in self.module_dict:
            self.module = module
        else:
            print(self.module_dict.keys())
            raise Exception(f"{module} is an Invalid Module Name Please Choose from List Above")

        module_name = self.module_dict[self.module]
        print(f' {module_name[3]} '.center(80, '*'))

        plot_title = module_name[3]
        xlim_range = [0.5, 5.5]
        marker_s = 40  # Default is 40
        alpha_val = 0.1  # Default is 0.1
        figure_dim = [10, 5]
        plot_colors = ['red', 'blue', 'green', 'orange', 'purple']
        
        plt.figure(figsize=(figure_dim[0], figure_dim[1]))

        combined_gsat = []
        combined_gmsl = []

        print(f'QUANTILES BY SSP')
        for scenario in range(len(self.scenarios)):
            this_scenario = f'{exp_name}{self.scenarios[scenario]}'

            # Pulls the FAIR Temperature GSAT data for the 19 year interval
            # Pulls the GMSL data out of a single module output
            gsat = self.get_gsat_data(f'{self.out_dir}/{this_scenario}.{self.module_dict["fair_gsat"][0]}')
            gmsl = self.get_module_data(filename=f'{self.out_dir}/{this_scenario}.{module_name[0]}')

            # Pulls the quantiles for all the gmsl ssps:
            quants = self.get_quants(gmsl, ssp_label=f'{this_scenario}')

            # Plots the data for the current SSP / Scenario
            plt.scatter(gsat, gmsl, marker='o', s=marker_s, color=plot_colors[scenario], alpha=alpha_val, edgecolors='none',
                        label=f'SSP{self.scenarios[scenario]}: {quants}')
            
            combined_gsat = np.append(combined_gsat, gsat)
            combined_gmsl = np.append(combined_gmsl, gmsl)
     
        # Overlays the binning data on the scatter plots
        ssp_comb = self.bin_data(combined_gsat,
                    combined_gmsl,
                    bin_start=0.25,
                    bin_stop=7,
                    interval=.5,
                    cutoff=200,
                    term_out=True,
                    plot_ax=False)


        plt.xlim(xlim_range[0], xlim_range[1])
        plt.ylim(module_name[1], module_name[2])


        plt.title(f'{plot_title} {self.version_flag}\n NSAMPS PER SCENARIO = {2000}')
        plt.xlabel('2081-2100 Average GSAT [C$^\circ$]')
        plt.ylabel(f'2100 GMSL [m]')
        plt.legend(bbox_to_anchor=(1,1), loc='upper left')
        plt.tight_layout()

        plt.savefig(f'{self.plot_dir}/{module}.png')
        #plt.show()