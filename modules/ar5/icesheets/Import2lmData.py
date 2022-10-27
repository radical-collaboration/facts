import numpy as np
import os
import sys
import h5py
import re
from netCDF4 import Dataset

'''
Import2lmData()

Imports two-layer model data set provided by Tim Hermans on 27 October 2020.  This
function reads in the netCDF file and filters the data for the variable, and scenario
of interest.  Additional filtering according to years, reference year(s), etc. will be
handled in another script.

'''
'''
def Import2lmData(variable="tmix", scenario="ssp585", directory="./"):

	# Open the netcdf file
	ncfiles = {"tmix": os.path.join(directory, "scmpy2LM_RCMIP_CMIP6calpm_ens_tmix.nc"), \
				"tdeep": os.path.join(directory, "scmpy2LM_RCMIP_CMIP6calpm_ens_tdeep.nc"), \
				"ohc": os.path.join(directory, "scmpy2LM_RCMIP_CMIP6calpm_ens_t_ohc.nc")}
	nc = Dataset(ncfiles[variable], "r")

	# Extract the years and scenarios
	years = nc.variables["years"][:]
	scenarios = nc.variables["scenario"][:]

	# Which scenario matches the one requested
	scenario_idx = np.flatnonzero(scenarios == scenario)

	# Extract the samples [samples, years]
	samps = np.squeeze(nc.variables[variable][scenario_idx,:,:])

	# Close the netcdf file
	nc.close()

	# Create the 2lm dictionary
	out_dict = {"samples": samps, "years": years}

	return(out_dict)


'''
'''
Import2lmData()

Imports two-layer model data set provided by Chris Smith on 30 October 2020.  This
function reads in the hdf5 file and filters the data for the variable, and scenario
of interest.  Additional filtering according to years, reference year(s), etc. will be
handled in another script.

'''
def Import2lmData(variable="surface_temperature", scenario="ssp585", directory="./", refyear_start=1995, refyear_end=2014, twinyear_start=2020, twinyear_end=2100,climate_fname="twolayer_SSPs.h5"):

	# Open the SSP hdf5 file
	sspfile = os.path.join(directory, climate_fname)
	hf = h5py.File(sspfile, 'r')

	# Do we have a temperature target scenario?
	scenario_test = re.search("^tlim(\d*\.?\d+)win(\d*\.?\d+)$", scenario)
	if(scenario_test):

		# Initialize samps array
		samps = []
		temp_samps = []

		# Loop through all available scenarios
		for this_scenario in hf.keys():

			# Skip "year" in this_scenario
			if this_scenario == "year":
				continue

			# Extract the samples
			try:
				these_samps = hf[this_scenario][variable][()]
				these_temp_samps = hf[this_scenario]["surface_temperature"][()]
			except:
				raise Exception("Cannot extract data for this combination: {} - {}".format(scenario, variable))

			# Append these samples to the main samps array
			samps.extend(these_samps.T)
			temp_samps.extend(these_temp_samps.T)

		# Make samps a numpy array
		samps = np.array(samps).T
		temp_samps = np.array(temp_samps).T

		# Get the years from the shape of the samps array
		years = np.arange(1750, 1750+samps.shape[0])

	# We have a standard SSP scenario
	else:

		try:
			samps = hf[scenario][variable][()]
		except:
			raise Exception("Cannot extract data for this combination: {} - {}".format(scenario, variable))
		years = np.arange(1750, 1750+samps.shape[0])

	# Calculate the reference period values (mean between refyear_start and refyear_end inclusive)
	year_idx = np.flatnonzero(np.logical_and(years >= refyear_start, years <= refyear_end))
	ref_vals = np.mean(samps[year_idx,:], axis=0)[None,:]
	samps = samps - ref_vals

	# If this is a temperature target scenario, filter for the samples within the target window
	# over the years 2020 - 2100
	if(scenario_test):

		# Extract the limit from the scenario string
		temp_target = float(scenario_test.group(1))
		temp_target_window = float(scenario_test.group(2))

		# Get the year indices for years 2020 - 2100
		window_years_idx = np.flatnonzero(np.logical_and(years >= twinyear_start, years <= twinyear_end))

		# Get the sample indices that match the filter
		samps_max = np.nanmax(temp_samps[window_years_idx,:], axis=0)
		match_idx = np.logical_and(samps_max >= temp_target-temp_target_window, samps_max <= temp_target+temp_target_window)

		# Subset samps for the samples that match the filter
		samps = samps[:,match_idx]

	# Close the input file
	hf.close()

	# Create the 2lm dictionary
	out_dict = {"samples": samps.T, "years": years}

	return(out_dict)




'''
Filter2lmData()

Filters a matrix of temperature trajectories produced by Import2lmData() for samples with
maximum temperature falling between two values over a span of time.

Returns a dictionary with "samples" = [matched_samples, years] and "years" = [years]

'''
def Filter2lmData(tlm_dict, filter_years=None, tmin=np.NINF, tmax=np.PINF):

	# Extract the data from the dict
	samps = tlm_dict["samples"]
	years = tlm_dict["years"]

	# Subset the data for the years requested
	if filter_years is not None:
		year_idx = np.isin(years, filter_years)
		samps = samps[:,year_idx]
		years = years[year_idx]

	# Find the maximum temperature for each sample
	samps_max = np.nanmax(samps, axis=1)

	# Which of these fall between the temperature limits
	match_idx = np.logical_and(samps_max >= tmin, samps_max <= tmax)

	# Fail if no samples match
	if not np.any(match_idx):
		raise Exception("Cannot find samples that match the target window. min = {}, max = {}".format(tmin, tmax))

	# Return the dictionary
	out_dict = {"samples": samps[match_idx,:], "years": years}
	return(out_dict)


if __name__ == "__main__":

	#varname = "ocean_heat_content"
	varname = "surface_temperature"
	scenario = "tlim1.5win0.25"
	#scenario = "ssp585"
	tlm_dict = Import2lmData(varname, scenario, "./", 1995, 2014)

	filtered_data_dict = Filter2lmData(tlm_dict, filter_years=np.arange(2000,2101))

	print(filtered_data_dict["samples"].shape)

	sys.exit()
