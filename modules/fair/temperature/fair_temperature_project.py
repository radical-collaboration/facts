import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
import fair
import argparse
import pickle
from datetime import datetime

from my_FAIR_forward import fair_scm as my_fair_scm


# Function that generates a sample from the FAIR model with provided emissions
def my_run_fair(args, emissions, reference_year=1750):
	"""
	# Run the rcmip emissions as-is, forward mode only.
	# This func never runs concentration-driven

	# emissions_dict should be like : {"ssp370": emissions array,
	#                                  "ssp245": emissions array,
	#                                  "ssp460": emissions array,
	#                                  "ssp585": emissions array }

	"""

	nt = 751  # 736

	#     with warnings.catch_warnings():
	#         warnings.simplefilter('ignore')
	F_solar = np.zeros(nt)
	F_volcanic = np.zeros(nt)
	natural = np.zeros((nt, 2))
	F_solar[:361] = args["F_solar"]
	F_volcanic[:361] = args["F_volcanic"]
	natural[:361, :] = args["natural"]
	natural[361:, :] = args["natural"][-1,:] # hold the last element constant for the rest of the array

	Cc, Fc, Tc, Tdeep, ariacic, lambda_effc, ohcc, heatfluxc = my_fair_scm(  #fair.forward.fair_scm(
			emissions=emissions,
			emissions_driven=True,
			C_pi=args["C_pi"].values, # new for round 2
			natural=natural,
			F_volcanic=F_volcanic,  # np.mean(F_volcanic),
			F_solar=F_solar,
			F2x=args["F2x"].values,
			r0=args["r0"].values,
			rt=args["rt"].values,
			rc=args["rc"].values,
			ariaci_out=True,
			ghg_forcing="Meinshausen",  # <- replace with args ?
			scale=np.asarray(args["scale"].values),
			aerosol_forcing="aerocom+ghan2",
			b_aero=np.asarray(args["b_aero"].values),
			ghan_params=np.asarray(args["ghan_params"].values),
			tropO3_forcing="thorhnill+skeie",  # <- replace with args ?
			ozone_feedback=args["ozone_feedback"].values,
			b_tro3=np.array(args["b_tro3"].values),
			E_pi=emissions[0, :],
			scaleAerosolAR5=False,
			scaleHistoricalAR5=False,
			fixPre1850RCP=False,
			aviNOx_frac=0,
			aCO2land=args["aCO2land"].values,
			stwv_from_ch4=args["stwv_from_ch4"].values,
			F_ref_BC=args["F_ref_BC"].values,
			E_ref_BC=args["E_ref_BC"].values,
			efficacy=np.ones(45),
			diagnostics="AR6",
			temperature_function="Geoffroy",
			lambda_global=args[
				"lambda_global"
			].values,  # this and the below only used in two-layer model
			deep_ocean_efficacy=args["deep_ocean_efficacy"].values,
			ocean_heat_capacity=np.asarray(args["ocean_heat_capacity"].values),
			ocean_heat_exchange=args["ocean_heat_exchange"].values,
		)

	# Return temperature and ocean heat content
	return(Tc, Tdeep, ohcc)


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



def fair_project_temperature(nsamps, seed, cyear_start, cyear_end, smooth_win, pipeline_id):

	# Load the preprocessed data
	preprocess_file = "{}_preprocess.pkl".format(pipeline_id)
	with open(preprocess_file, 'rb') as f:
		preprocess_data = pickle.load(f)

	emis = preprocess_data["emis"]
	REFERENCE_YEAR = preprocess_data["REFERENCE_YEAR"]
	scenario = preprocess_data["scenario"]
	rcmip_file = preprocess_data["rcmip_file"]

	# Load the fit data
	fit_file = "{}_fit.pkl".format(pipeline_id)
	with open(fit_file, 'rb') as f:
		fit_data = pickle.load(f)

	pars = fit_data["pars"]
	param_file = fit_data["param_file"]

	# Initialize lists to hold samples
	temps = []
	deeptemps = []
	ohcs = []

	# Projection years
	proj_years = np.arange(REFERENCE_YEAR, 2501)

	# What version of FAIR are we using?
	fair_version = fair.__version__

	# How many parameter sets do we have?
	nsims = len(pars["simulation"])

	# Generate nsamps of simulation indices to sample
	rng = np.random.default_rng(seed)
	if nsamps > nsims:
		run_idx = np.arange(nsims)
		sample_idx = rng.choice(nsims, nsamps, nsamps>nsims)
	else:
		run_idx = rng.choice(nsims, nsamps, nsamps>nsims)
		sample_idx = np.arange(nsamps)

	# Run the FAIR model
	for i in run_idx:
		this_pars = pars.isel(simulation=i)
		this_temp, this_deeptemp, this_ohc = my_run_fair(this_pars, emis)
		temps.append(this_temp)
		deeptemps.append(this_deeptemp)
		ohcs.append(this_ohc)

	# Recast the output as numpy arrays
	temps = np.array(temps)
	deeptemps = np.array(deeptemps)
	ohcs = np.array(ohcs)

	# Center and smooth the samples
	temps = CenterSmooth(temps, proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
	deeptemps = CenterSmooth(deeptemps, proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)
	ohcs = CenterSmooth(ohcs, proj_years, cyear_start=cyear_start, cyear_end=cyear_end, smooth_win=smooth_win)

	# Conform the output to shapes appropriate for output
	temps = temps[sample_idx,:,np.newaxis]
	deeptemps = deeptemps[sample_idx,:,np.newaxis]
	ohcs = ohcs[sample_idx,:,np.newaxis]

	# Set up the global attributes for the output netcdf files
	attrs = {"Source": "FACTS",
			 "Date Created": str(datetime.now()),
			 "Description": (
				 f"Fair v={fair_version} scenario simulations with AR6-calibrated settings."
			 " Simulations based on parameters developed here: https://github.com/chrisroadmap/ar6/tree/main/notebooks."
			 " Parameters obtained from: https://zenodo.org/record/5513022#.YVW1HZpByUk."),
			"Method": (
				"Temperature and ocean heat content were returned from fair.foward.fair_scm() in emission-driven mode."),
			"Scenario emissions file": rcmip_file,
			"FAIR Parameters file": param_file,
			"FaIR version": fair_version,
			 "Scenario": scenario,
			 "Centered": "{}-{} mean".format(cyear_start, cyear_end),
			 "Smoothing window": "{} years".format(smooth_win),
			 "Note": "Code provided by Kelly McCusker of Rhodium Group Climate Impact Lab and adapted for use in FACTS."
			}

	# Create the variable datasets
	tempds = xr.Dataset({"surface_temperature": (("samples", "years", "locations"), temps, {"units":"degC"}),
							"lat": (("locations"), [np.inf]),
							"lon": (("locations"), [np.inf])},
		coords={"years": proj_years, "locations": [-1], "samples": np.arange(nsamps)}, attrs=attrs)

	deeptempds = xr.Dataset({"deep_ocean_temperature": (("samples", "years", "locations"), deeptemps, {"units":"degC"}),
							"lat": (("locations"), [np.inf]),
							"lon": (("locations"), [np.inf])},
		coords={"years": proj_years, "locations": [-1], "samples": np.arange(nsamps)}, attrs=attrs)

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
							"deep_ocean_temperature": (("years","samples"), deeptemps[::,::,0].transpose(), {"units":"degC"}),
							"ocean_heat_content": (("years","samples"), ohcs[::,::,0].transpose(), {"units":"J"})},
		coords={"years": proj_years, "samples": np.arange(nsamps)}, attrs=attrs)
	pooledds.to_netcdf("{}_climate.nc".format(pipeline_id), group=scenario,encoding={"ocean_heat_content": {"dtype": "float32", "zlib": True, "complevel":4},
		"surface_temperature": {"dtype": "float32", "zlib": True, "complevel":4},
		"deep_ocean_temperature": {"dtype": "float32", "zlib": True, "complevel":4}})
	yearsds = xr.Dataset({"year": proj_years})
	yearsds.to_netcdf("{}_climate.nc".format(pipeline_id), mode='a')

	# Done
	return(None)


if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the project stage for the FAIR AR6 temperature module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--nsamps', help="Number of samples to create (uses replacement if nsamps > n parameters) (default=10)", type=int, default=10)
	parser.add_argument('--seed', help="Seed value for random number generator (default=1234)", type=int, default=1234)
	parser.add_argument('--cyear_start', help="Start year of temporal range for centering (default=1850)", type=int, default=1850)
	parser.add_argument('--cyear_end', help="End year of temporal range for centering (default=1900)", type=int, default=1900)
	parser.add_argument('--smooth_win', help="Number of years to use as a smoothing window (default=19)", type=int, default=19)

	# Parse the arguments
	args = parser.parse_args()

	# Run the code
	fair_project_temperature(nsamps=args.nsamps, seed=args.seed, cyear_start=args.cyear_start, cyear_end=args.cyear_end, smooth_win=args.smooth_win, pipeline_id=args.pipeline_id)


	sys.exit()
