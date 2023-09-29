import argparse
import os
import sys
import pickle as p
import numpy as np
import scipy.io
import h5py

''' bamber19_preprocess_icesheets.py

This runs the preprocessing stage for the Bamber et al. 2019 ice sheet component of the IPCC AR6
workflow.

Parameters:
pipeline_id = Unique identifier for the pipeline running this code

'''

def bamber19_preprocess_icesheets(pyear_start, pyear_end, pyear_step, baseyear, scenario, pipeline_id):
	
	scenario_map = {"rcp85": 'corefileH', "rcp26": 'corefileL',  \
					"tlim2.0win0.25": 'corefileL', "tlim5.0win0.25": 'corefileH'}
	this_corefile = scenario_map[scenario]
	
	# Define the target years
	targyears = np.arange(pyear_start, pyear_end+1, pyear_step)
	
	# Load the data
	filename = os.path.join(os.path.dirname(__file__), "SLRProjections190726core_SEJ_full.mat")
	mat = scipy.io.loadmat(filename)
	
	wais_samps, eais_samps, gis_samps = ExtractSamples(mat, this_corefile, targyears, baseyear)
	OutputData(pipeline_id, eais_samps, wais_samps, gis_samps, scenario, targyears, baseyear)

def bamber19_preprocess_icesheets_temperaturedriven(pyear_start, pyear_end, pyear_step, baseyear, scenario, climate_data_file, pipeline_id):
	# Load the temperature data
	SAT,Time,NumTensemble = GetSATData(climate_data_file, scenario)

	# Rest of code needs to be written
	print('This does not work yet')

	wais_sampsH, eais_sampsH, gis_sampsH = ExtractSamples(mat, "corefileH", targyears, baseyear)
	wais_sampsL, eais_sampsL, gis_sampsL = ExtractSamples(mat, "corefileL", targyears, baseyear)

	# For each temperature time series sample, compute a score in 2050 and 2100 that compares to the Bamber
	# low and high scenarios (in 2050, low=1.5 C and high=2.0 C; in 2100, 2 C and 5 C, all relative to 1850-1900)
	# -- ie so score = 0 if below Low and 1 of above High, and linearly interpolated in between.
	# Average together 2050 and 2100 scores, then use score as probability of drawing from the High distribution
	# (as opposed to Low distribution) for that time series.

	# for each temperature sample you pick which scenario to use based on comparing the weighting factor to a random draw

	OutputData(pipeline_id, eais_samps, wais_samps, gis_samps, scenario, targyears, baseyear)

def GetSATData(fname, scenario, refyear_start=1850, refyear_end=1900, year_start=1900, year_end=2300):

	# Open the file
	df_ssp = h5py.File(fname,'r')

	# Extract the surface temperature for this scenario
	if scenario not in df_ssp.keys():
		raise Except("Scenario {} not found in {}".format(scenario, fname))
	ssp_folder = df_ssp.get(scenario)
	sat_ssp = ssp_folder.get('surface_temperature')

	# Get the number of ensemble members from the data
	_,nens = sat_ssp.shape

	# Extract the years available
	sat_years = df_ssp['year'][()]

	# Which indices align with the reference and trim years
	refyear_start_idx = np.flatnonzero(sat_years == refyear_start)[0]
	refyear_end_idx = np.flatnonzero(sat_years == refyear_end)[0] #+ 1
	year_start_idx = np.flatnonzero(sat_years == year_start)[0]
	year_end_idx = np.flatnonzero(sat_years == year_end)[0] + 1

	# Normalize and crop temperature series
	Time = np.arange(year_start, year_end+1)
	SATave = np.mean(sat_ssp[refyear_start_idx:refyear_end_idx,:], axis=0)
	SAT = sat_ssp[year_start_idx:year_end_idx,:] - SATave

	# Close the h5 file
	df_ssp.close()

	# Done
	return(SAT, Time, nens)

def ExtractSamples(mat, this_corefile, targyears, baseyear):
	# Get the years available from the matlab core file
	mat_years = np.squeeze(mat[this_corefile][0,0][27])
	
	# Determine which matlab year indices match our target years
	mat_years_idx = np.isin(mat_years, targyears)
	
	# Get the samples from the matlab core file
	samps = mat[this_corefile][0,0][21]
	
	# Extract the samples for each ice sheet
	eais_samps = samps[:,20,:]
	wais_samps = samps[:,19,:]
	gis_samps = samps[:,18,:]
	
	# Get the values for the baseyear of interest
	eais_refs = np.apply_along_axis(FindRefVals, axis=1, arr=eais_samps, years=mat_years, baseyear=baseyear)
	wais_refs = np.apply_along_axis(FindRefVals, axis=1, arr=wais_samps, years=mat_years, baseyear=baseyear)
	gis_refs = np.apply_along_axis(FindRefVals, axis=1, arr=gis_samps, years=mat_years, baseyear=baseyear)
	
	# Center the projections to the reference period
	eais_samps -= eais_refs[:,np.newaxis]
	wais_samps -= wais_refs[:,np.newaxis]
	gis_samps -= gis_refs[:,np.newaxis]
	
	# Subset for the target years
	eais_samps = eais_samps[:,mat_years_idx]
	wais_samps = wais_samps[:,mat_years_idx]
	gis_samps = gis_samps[:,mat_years_idx]

	return wais_samps, eais_samps, gis_samps


def OutputData(pipeline_id, eais_samps, wais_samps, gis_samps, scenario, targyears, baseyear):

	# Sum up the components to get total AIS samples
	ais_samps = eais_samps + wais_samps

	# Populate the output dictionary
	outdata = {'eais_samps': eais_samps, 'wais_samps': wais_samps, 'ais_samps': ais_samps, \
				'gis_samps': gis_samps, 'scenario': scenario, 'targyears': targyears, 'baseyear': baseyear}
	
	# Define the data directory
	outdir = os.path.dirname(__file__)

	# Write the rates data to a pickle file
	outfile = open(os.path.join(outdir, "{}_data.pkl".format(pipeline_id)), 'wb')
	p.dump(outdata, outfile)
	outfile.close()


def FindRefVals(timeseries, years, baseyear):
	
	# Append a zero to the beginning of the timeseries at year 2000
	timeseries = np.append(np.array([0.0]), timeseries)
	years = np.append(np.array([2000]), years)
	
	# Interpolate to the appropriate base year
	ref_val = np.interp(baseyear, years, timeseries, left=0.0)
	
	# Return the value
	return(ref_val)


if __name__ == '__main__':	
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the pre-processing stage for the IPCC AR6 Bamber et al. 2019 ice sheet workflow",\
	epilog="Note: This is meant to be run as part of the IPCC AR6 module within the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2100, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--scenario', help="Emissions scenario of interest [default=rcp85]", choices=["rcp85", "rcp26", "tlim2.0win0.25", "tlim5.0win0.25"], default="rcp85")
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str, default="")
	parser.add_argument('--baseyear', help="Year to which projections are referenced [default = 2000]", default=2000, type=int)
	
	# Parse the arguments
	args = parser.parse_args()
	
	if len(args.climate_data_file) == 0:
		# Run the preprocessing stage with the user defined RCP scenario
		bamber19_preprocess_icesheets(args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario, args.pipeline_id)
	else:
		bamber19_preprocess_icesheets_temperaturedriven(args.pyear_start, args.pyear_end, args.pyear_step, args.baseyear, args.scenario, args.climate_data_file, args.pipeline_id)

	# Done
	sys.exit()