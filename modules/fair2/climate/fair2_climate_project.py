import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr
from fill_from_rcmip import fill_from_rcmip
from fill_from import *

# Function that prepares the RCMIP emissions data
def prep_rcmip_emissions_conc(scenario, rcmip_file):
	pass


def fair2_project_climate(scenario,rcmip_file, concentration_file,forcing_file, calibration_file, species_config_file, volcanic_erf, solar_erf, reference_year, pipeline_id, nsamps,cyear_start, cyear_end,
						  pyear_start,pyear_end,smooth_win,seed):
	
	from fair import FAIR
	from fair.interface import fill, initialise
	from fair.io import read_properties
	
	f = FAIR(ch4_method="Thornhill2021")

	# This should not be hardcoded
	f.define_time(reference_year, 2300, 1)  # start, end, step

	scenarios = [scenario]
	f.define_scenarios(scenarios)

	df_configs = pd.read_csv(calibration_file, index_col=0)
	configs = df_configs.index  # this is used as a label for the "config" axis
	f.define_configs(configs)
	configs
	df_configs.head()

	species, properties = read_properties(filename=species_config_file)
	f.define_species(species, properties)
	
	f.allocate()

	#f.fill_from_rcmip()
	if rcmip_file=='./rcmip/rcmip-emissions-annual-means-v5-1-0.csv':
		fill_from_rcmip(f, rcmip_file=rcmip_file)
	else:
		fill_from_csv(f,emissions_file=rcmip_file, concentration_file=concentration_file, forcing_file=forcing_file)
	f.emissions

	df_emis = pd.read_csv(rcmip_file)
	gfed_sectors = [
    	"Emissions|NOx|MAGICC AFOLU|Agricultural Waste Burning",
    	"Emissions|NOx|MAGICC AFOLU|Forest Burning",
    	"Emissions|NOx|MAGICC AFOLU|Grassland Burning",
    	"Emissions|NOx|MAGICC AFOLU|Peat Burning",
	]
	for scenario in scenarios:
		f.emissions.loc[dict(specie="NOx", scenario=scenario)] = (
        	df_emis.loc[
            (df_emis["Scenario"] == scenario)
            & (df_emis["Region"] == "World")
            & (df_emis["Variable"].isin(gfed_sectors)),
            "1750":"2300",
        ]
        .interpolate(axis=1)
        .values.squeeze()
        .sum(axis=0)
        * 46.006
        / 30.006
        + df_emis.loc[
            (df_emis["Scenario"] == scenario)
            & (df_emis["Region"] == "World")
            & (df_emis["Variable"] == "Emissions|NOx|MAGICC AFOLU|Agriculture"),
            "1750":"2300",
        ]
        .interpolate(axis=1)
        .values.squeeze()
        + df_emis.loc[
            (df_emis["Scenario"] == scenario)
            & (df_emis["Region"] == "World")
            & (df_emis["Variable"] == "Emissions|NOx|MAGICC Fossil and Industrial"),
            "1750":"2300",
        ]
        .interpolate(axis=1)
        .values.squeeze()
    )[:550, None]

	df_solar = pd.read_csv(solar_erf, index_col="year")
	df_volcanic = pd.read_csv(volcanic_erf)

	solar_forcing = np.zeros(551)
	volcanic_forcing = np.zeros(551)
	volcanic_forcing[:352] = df_volcanic.erf.values
	solar_forcing = df_solar["erf"].loc[reference_year:2300].values

	trend_shape = np.ones(551)
	trend_shape[:271] = np.linspace(0, 1, 271)

	fill(
    f.forcing,
    volcanic_forcing[:, None, None] * df_configs["fscale_Volcanic"].values.squeeze(),
    specie="Volcanic",)
	
	fill(
    f.forcing,
    solar_forcing[:, None, None] * df_configs["fscale_solar_amplitude"].values.squeeze()
    + trend_shape[:, None, None] * df_configs["fscale_solar_trend"].values.squeeze(),
    specie="Solar",)

	fill(f.climate_configs["ocean_heat_capacity"], df_configs.loc[:, "clim_c1":"clim_c3"].values)
	fill(f.climate_configs["ocean_heat_transfer"],df_configs.loc[:, "clim_kappa1":"clim_kappa3"].values,)
	fill(f.climate_configs["deep_ocean_efficacy"], df_configs["clim_epsilon"].values.squeeze())
	fill(f.climate_configs["gamma_autocorrelation"], df_configs["clim_gamma"].values.squeeze())
	fill(f.climate_configs["sigma_eta"], df_configs["clim_sigma_eta"].values.squeeze())
	fill(f.climate_configs["sigma_xi"], df_configs["clim_sigma_xi"].values.squeeze())
	fill(f.climate_configs["seed"], df_configs["seed"])
	fill(f.climate_configs["stochastic_run"], True)
	fill(f.climate_configs["use_seed"], True)
	fill(f.climate_configs["forcing_4co2"], df_configs["clim_F_4xCO2"])

	f.fill_species_configs(filename=species_config_file)

	# carbon cycle
	fill(f.species_configs["iirf_0"], df_configs["cc_r0"].values.squeeze(), specie="CO2")
	fill(f.species_configs["iirf_airborne"], df_configs["cc_rA"].values.squeeze(), specie="CO2")
	fill(f.species_configs["iirf_uptake"], df_configs["cc_rU"].values.squeeze(), specie="CO2")
	fill(f.species_configs["iirf_temperature"], df_configs["cc_rT"].values.squeeze(), specie="CO2")

	# aerosol indirect
	fill(f.species_configs["aci_scale"], df_configs["aci_beta"].values.squeeze())
	fill(f.species_configs["aci_shape"], df_configs["aci_shape_so2"].values.squeeze(), specie="Sulfur")
	fill(f.species_configs["aci_shape"], df_configs["aci_shape_bc"].values.squeeze(), specie="BC")
	fill(f.species_configs["aci_shape"], df_configs["aci_shape_oc"].values.squeeze(), specie="OC")

	# aerosol direct
	for specie in [
    	"BC", 
    	"CH4", 
    	"N2O",
    	"NH3", 
    	"NOx",
    	"OC", 
    	"Sulfur", 
    	"VOC",
    	"Equivalent effective stratospheric chlorine"
	]:
		fill(f.species_configs["erfari_radiative_efficiency"], df_configs[f"ari_{specie}"], specie=specie)

	# forcing scaling
	for specie in [
    	"CO2", 
    	"CH4", 
    	"N2O", 
    	"Stratospheric water vapour",
    	"Contrails", 
    	"Light absorbing particles on snow and ice", 
    	"Land use"
	]:
		fill(f.species_configs["forcing_scale"], df_configs[f"fscale_{specie}"].values.squeeze(), specie=specie)

	# the halogenated gases all take the same scale factor
	for specie in [
    	"CFC-11",
    	"CFC-12",
    	"CFC-113",
    	"CFC-114",
    	"CFC-115",
    	"HCFC-22",
    	"HCFC-141b",
    	"HCFC-142b",
    	"CCl4",
    	"CHCl3",
    	"CH2Cl2",
    	"CH3Cl",
    	"CH3CCl3",
    	"CH3Br",
    	"Halon-1211",
    	"Halon-1301",
    	"Halon-2402",
    	"CF4",
    	"C2F6",
    	"C3F8",
    	"c-C4F8",
    	"C4F10",
    	"C5F12",
    	"C6F14",
    	"C7F16",
    	"C8F18",
    	"NF3",
    	"SF6",
    	"SO2F2",
    	"HFC-125",
    	"HFC-134a",
    	"HFC-143a",
    	"HFC-152a",
    	"HFC-227ea",
    	"HFC-23",
    	"HFC-236fa",
    	"HFC-245fa",
    	"HFC-32",
    	"HFC-365mfc",
    	"HFC-4310mee",
		]:
		fill(f.species_configs["forcing_scale"], df_configs["fscale_minorGHG"].values.squeeze(), specie=specie)

	# ozone
	for specie in ["CH4", "N2O", "Equivalent effective stratospheric chlorine", "CO", "VOC", "NOx"]:
		fill(f.species_configs["ozone_radiative_efficiency"], df_configs[f"o3_{specie}"], specie=specie)

	# initial value of CO2 concentration (but not baseline for forcing calculations)
	fill(f.species_configs["baseline_concentration"], df_configs["cc_co2_concentration_1750"].values.squeeze(), specie="CO2")

	initialise(f.concentration, f.species_configs["baseline_concentration"])
	initialise(f.forcing, 0)
	initialise(f.temperature, 0)
	initialise(f.cumulative_emissions, 0)
	initialise(f.airborne_emissions, 0)

	f.run()

	ohc = f.ocean_heat_content_change
	temp3 = f.temperature.loc[dict(layer=[0, 1, 2])]

	rng = np.random.default_rng(seed)
	nsims = 1001
	si = 0

	if nsamps > nsims:
		run_idx = np.arange(nsims)
		sample_idx = rng.choice(nsims, nsamps, nsamps>nsims)
	else:
		run_idx = rng.choice(nsims, nsamps, nsamps>nsims)
		sample_idx = np.arange(nsamps)

	proj_years = np.arange(reference_year,2301)
	
	
	temp3t = np.transpose(np.array(temp3[:,0,sample_idx,0]))
	ohct = np.transpose(np.array(ohc[:,0,sample_idx]))

    
	temps = np.empty((nsamps, temp3t.shape[1], 1))
	deeptemps = np.empty((nsamps, temp3t.shape[1], 2, 1))
	deeptempst =  np.empty((temp3t.shape[1], nsamps,  2, 1))
	ohcs = np.empty((nsamps, temp3t.shape[1], 1))
		
	temps[:,:,0] = temp3t
	temps = np.array(temps)
	ohcs[:,:,0] = np.array(ohct)
	temps = np.array(temps)
	deeptemps[:,:,0,0] = np.transpose(np.array(temp3[:,0,sample_idx,1]))
	deeptemps[:,:,1,0] = np.transpose(np.array(temp3[:,0,sample_idx,2]))
    
   	# deeptemps = np.array(deeptemps)
   	# ohcs = np.array(ohcs)
   
   	# Center and smooth the samples
	temps = CenterSmooth(temps, proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
	deeptemps[:,:,0,0] = CenterSmooth(deeptemps[:,:,0,0], proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
	deeptemps[:,:,1,0] = CenterSmooth(deeptemps[:,:,1,0], proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
	deeptempst[:,:,0,0] = np.transpose(deeptemps[:,:,0,0])
	deeptempst[:,:,1,0] = np.transpose(deeptemps[:,:,1,0])
		
	ohcs = CenterSmooth(ohcs, proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
   
   	# Conform the output to shapes appropriate for output
   	# deeptemps = deeptemps[sample_idx,:,np.newaxis]
   	# ohcs = ohcs[sample_idx,:,np.newaxis]
   
   	# Set up the global attributes for the output netcdf files
	attrs = {"Source": "FACTS",
   			 "Date Created": str(datetime.now()),
   			 "Description": (
   				 "Fair v=2.1.2 scenario simulations with AR6-calibrated settings."
   			 " Simulations based on parameters developed here: https://github.com/chrisroadmap/ar6/tree/main/notebooks."
   			 " Parameters obtained from: https://zenodo.org/records/8399112."),
   			"Method": (
   				"Temperature retunr from f = FAIR(ch4_method='Thornhill2021')."),
   			"Scenario emissions file": rcmip_file,
   			"FAIR Parameters file": calibration_file,
   			"FaIR version": '2.1.6',
   			 "Scenario": scenario,
   			 "Centered": "{}-{} mean".format(cyear_start, cyear_end),
   			 "Smoothing window": "{} years".format(smooth_win),
   			 "Note": "Code Provided by Víctor Malagón Santos and adapted to FACTS by Alex Reedy."
   			}
   
		# Create the variable datasets
	tempds = xr.Dataset({"surface_temperature": (("samples", "years", "locations"), temps, {"units":"degC"}),
   							"lat": (("locations"), [np.inf]),
   							"lon": (("locations"), [np.inf])},
   		coords={"years": proj_years, "locations": [-1], "samples": np.arange(nsamps)}, attrs=attrs)
   
	deeptempds = xr.Dataset({"deep_ocean_temperature": (("samples", "years", "layers", "locations"), deeptemps, {"units":"degC"}),
   							"lat": (("locations"), [np.inf]),
   							"lon": (("locations"), [np.inf])},
   	 	coords={"years": proj_years, "locations": [-1], "samples": np.arange(nsamps), "layers": [1,2]}, attrs=attrs)
   
	ohcds = xr.Dataset({"ocean_heat_content": (("samples", "years", "locations"), ohcs, {"units":"J"}),
   							"lat": (("locations"), [np.inf]),
   							"lon": (("locations"), [np.inf])},
   	 	coords={"years": proj_years, "locations": [-1], "samples": np.arange(nsamps)}, attrs=attrs)
   
		# Write the datasets to netCDF
	tempds.to_netcdf("{}_gsat.nc".format(pipeline_id), encoding={"surface_temperature": {"dtype": "float32", "zlib": True, "complevel":4}})
	deeptempds.to_netcdf("{}_oceantemp.nc".format(pipeline_id), encoding={"deep_ocean_temperature": {"dtype": "float32", "zlib": True, "complevel":4}})
	ohcds.to_netcdf("{}_ohc.nc".format(pipeline_id), encoding={"ocean_heat_content": {"dtype": "float32", "zlib": True, "complevel":4}})
   
		# create a single netCDF file that is compatible with modules expecting parameters organized in a certain fashion
	pooledds = xr.Dataset({"surface_temperature": (("years","samples"), temps[::,::,0].transpose(), {"units":"degC"}),
   							"deep_ocean_temperature": (("years","samples", "layers"), deeptempst[::,::,::,0], {"units":"degC"}),
   							"ocean_heat_content": (("years","samples"), ohcs[::,::,0].transpose(), {"units":"J"})},
   	 	coords={"years": proj_years, "samples": np.arange(nsamps), "layers": [1,2]}, attrs=attrs)

	pooledds.to_netcdf("{}_climate.nc".format(pipeline_id), group=scenario,encoding={"ocean_heat_content": {"dtype": "float32", "zlib": True, "complevel":4},
   	 	"surface_temperature": {"dtype": "float32", "zlib": True, "complevel":4},
   	 	"deep_ocean_temperature": {"dtype": "float32", "zlib": True, "complevel":4}})
	yearsds = xr.Dataset({"year": proj_years})
    
	yearsds.to_netcdf("{}_climate.nc".format(pipeline_id), mode='a')
	
	return(None)

# Function that centers and smooths a set of temperature trajectories
def CenterSmooth(temp, years, cyear_start, cyear_end, smooth_win):

	# Which years fall into the range of years for centering?
	ref_idx = np.flatnonzero(np.logical_and(years >= cyear_start, years <= cyear_end))

	# Take the mean over the centering range for each sample and subtract
	ref_temp = np.mean(temp[:,ref_idx], axis=1)
	center_temp = temp - ref_temp[:,np.newaxis]

	# Smooth the centered temperatures
	if smooth_win > 1:
		csmooth_temp = np.apply_along_axis(Smooth, axis=1, arr=center_temp, w=smooth_win)
	else:
		csmooth_temp = center_temp

	# Return the centered and smoothed temperature trajectories
	return(csmooth_temp)


# Function that smooths data in similar fashion to Matlab's "smooth" function
def Smooth(x, w=5):
	out0 = np.convolve(x, np.ones(w,dtype='double'), 'valid')/w
	r = np.arange(1,w-1,2, dtype="double")
	start = np.cumsum(x[:w-1])[::2]/r
	stop = (np.cumsum(x[:-w:-1])[::2]/r)[::-1]
	y = np.concatenate((start, out0, stop))
	return(y)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the preprocessing stage for the FAIR AR6 temperature module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--rcmip_file', help="Full path to RCMIP emissions file", default="./rcmip/rcmip-emissions-annual-means-v5-1-0.csv")
	parser.add_argument('--concentration_file', help='Full path to concentration file', default="./rcmip/rcmip-concentrations-annual-means-v5-1-0.csv")
	parser.add_argument('--forcing_file', help='Full path to forcing file', default="./rcmip/rcmip-radiative-forcing-annual-means-v5-1-0.csv")
	parser.add_argument('--calibration_file', help="Full path to the climate calibration params file", default="./parameters/calibrated_constrained_parameters.csv")
	parser.add_argument('--species_config_file', help="Full path to the species configuration file", default="./parameters/species_configs_properties_calibration1.2.0.csv")
	parser.add_argument('--volcanic_erf', help='Full path to the Volcanic ERF timebounds file', default="./parameters/volcanic_ERF_1750-2101_timebounds.csv")
	parser.add_argument('--solar_erf', help='Full path to the Solar ERF timebounds file', default="./parameters/solar_erf_timebounds.csv")
	parser.add_argument('--scenario', help="SSP Emissions scenario",  default="ssp585")
	parser.add_argument('--reference_year', help='Reference Year', default=1750)
	parser.add_argument('--nsamps', help="Number of samples to create (uses replacement if nsamps > n parameters) (default=10)", type=int, default=10)
	parser.add_argument('--seed', help="Seed value for random number generator (default=1234)", type=int, default=1234)
	parser.add_argument('--cyear_start', help="Start year of temporal range for centering (default=1850)", type=int, default=1850)
	parser.add_argument('--cyear_end', help="End year of temporal range for centering (default=1900)", type=int, default=1900)
	parser.add_argument('--smooth_win', help="Number of years to use as a smoothing window (default=19)", type=int, default=19)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	fair2_project_climate(
		scenario=args.scenario, 
		rcmip_file=args.rcmip_file, 
		concentration_file=args.concentration_file,
		forcing_file=args.forcing_file,
		calibration_file=args.calibration_file, 
		species_config_file=args.species_config_file,
		volcanic_erf=args.volcanic_erf,
		solar_erf=args.solar_erf,
		pipeline_id=args.pipeline_id,
		reference_year=args.reference_year,
		nsamps=args.nsamps,
		seed=args.seed,
		cyear_start=args.cyear_start,
		cyear_end=args.cyear_end,
		smooth_win=args.smooth_win,
		pyear_start=args.pyear_start,
		pyear_end=args.pyear_end
		)

	sys.exit()