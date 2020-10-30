import numpy as np
import os
import sys
from netCDF4 import Dataset

'''
Import2lmData()

Imports two-layer model data set provided by Tim Hermans on 27 October 2020.  This 
function reads in the netCDF file and filters the data for the variable, and scenario 
of interest.  Additional filtering according to years, reference year(s), etc. will be 
handled in another script.

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
	
	varname = "tmix"
	scenario = "rcp85"
	tlm_dict = Import2lmData(varname, scenario)
	
	filtered_data_dict = Filter2lmData(tlm_dict, filter_years=np.arange(2000,2101), tmin=4.5, tmax=5.5)
	
	print(np.unique(filtered_data_dict["samples"], axis=0))
	
	sys.exit()