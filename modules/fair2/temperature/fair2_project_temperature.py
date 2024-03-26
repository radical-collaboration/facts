import os
import sys
import argparse
import pickle
import numpy as np
import pandas as pd
import pooch


# Function that prepares the RCMIP emissions data
def prep_rcmip_emissions_conc(scenario, rcmip_file):
	pass


def fair2_preprocess_temperature(scenario,rcmip_file, pipeline_id):
	from fair import FAIR
	from fair.interface import fill, initialise
	from fair.io import read_properties
	
	f = FAIR(ch4_method="Thornhill2021")

	f.define_time(1750, 2300, 1)  # start, end, step

	scenarios = [scenario]
	f.define_scenarios(scenarios)

	fair_params_1_2_0_obj = pooch.retrieve(
    url = 'https://zenodo.org/record/8399112/files/calibrated_constrained_parameters.csv',
    known_hash = 'md5:de3b83432b9d071efdd1427ad31e9076')
	df_configs = pd.read_csv(fair_params_1_2_0_obj, index_col=0)
	configs = df_configs.index  # this is used as a label for the "config" axis
	f.define_configs(configs)
	configs
	df_configs.head()

	species, properties = read_properties(filename='species_configs_properties_calibration1.2.0.csv')
	f.define_species(species, properties)
	
	f.allocate()

	f.fill_from_rcmip()
	f.emissions

	rcmip_emissions_file = pooch.retrieve(
    url="doi:10.5281/zenodo.4589756/rcmip-emissions-annual-means-v5-1-0.csv",
    known_hash="md5:4044106f55ca65b094670e7577eaf9b3",
)
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
		
	solar_obj = pooch.retrieve(
    url = 'https://raw.githubusercontent.com/chrisroadmap/fair-add-hfc/main/data/solar_erf_timebounds.csv',
    known_hash = 'md5:98f6f4c5309d848fea89803683441acf',
	)

	volcanic_obj = pooch.retrieve(
    url = 'https://raw.githubusercontent.com/chrisroadmap/fair-calibrate/main/data/forcing/volcanic_ERF_1750-2101_timebounds.csv',
    known_hash = 'md5:c0801f80f70195eb9567dbd70359219d',
	)	

	df_solar = pd.read_csv(solar_obj, index_col="year")
	df_volcanic = pd.read_csv(volcanic_obj)

	solar_forcing = np.zeros(551)
	volcanic_forcing = np.zeros(551)
	volcanic_forcing[:352] = df_volcanic.erf.values
	solar_forcing = df_solar["erf"].loc[1750:2300].values

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

	f.fill_species_configs(filename='species_configs_properties_calibration1.2.0.csv')

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
	
	return(None)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the preprocessing stage for the FAIR AR6 temperature module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--rcmip_file', help="Full path to RCMIP emissions file", default="./rcmip/rcmip-emissions-annual-means-v5-1-0.csv")
	parser.add_argument('--scenario', help="SSP Emissions scenario",  default="ssp585")

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	fair2_preprocess_temperature(scenario=args.scenario, rcmip_file=args.rcmip_file, pipeline_id=args.pipeline_id)


	sys.exit()