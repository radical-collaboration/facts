import numpy as np
import argparse
import pickle
import os
import re
import time
import sys
from netCDF4 import Dataset
import warnings

''' ipccar6_project_larmipicesheet.py

This is the projection stage for the LARMIP ice sheet component of the IPCC AR6 module set.

Parameters:
nsamps              Number of samples to produce
replace             Allow sampling with replacement
rngseed             Seed for the random number generator
pipeline_id         Unique identifier to attach to this pipeline


Output:
"{pipeline_id}_projections.pkl" = User requested samples
"{pipeline_id}_(AIS|EAIS|WAIS|PEN)_globalsl.nc" = Sampled global projections in netCDF file

'''

def ipccar6_project_larmipicesheet(nsamps, pipeline_id, cyear_start, cyear_end, replace, rngseed):

	# Load the data file
	datafilename = "{}_data.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	years = my_data["targyears"]
	scenario = my_data["scenario"]
	baseyear = my_data["baseyear"]
	ais_samples = my_data["ais_samps"]
	eais_samples = my_data["eais_samps"]
	wais_samples = my_data["wais_samps"]
	pen_samples = my_data["pen_samps"]

	# Load the smb preprocessed data
	datafilename = "{}_smbdata.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_data = pickle.load(f)

	temp_mean = my_data['temp_mean']
	temp_sd = my_data['temp_sd']
	inttemp_mean = my_data['inttemp_mean']
	inttemp_sd = my_data['inttemp_sd']
	#temp_samples = my_data['temp_samples']
	data_years = my_data['data_years']

	# Load the fit data
	datafilename = "{}_fit.pkl".format(pipeline_id)
	datafile = os.path.join(os.path.dirname(__file__), datafilename)

	with open(datafile, 'rb') as f:
		my_fit = pickle.load(f)

	# How many larmip samples do we take from each model?
	(nmodels, nsamples, nyears) = eais_samples.shape
	nsamps_per_model = np.repeat(int(nsamps / nmodels), nmodels)

	# Divy up the remaining larmip samples across randomly chosen models
	rng = np.random.default_rng(rngseed)
	nsamps_remain = nsamps % nmodels
	rand_model_idx = rng.choice(nmodels, size=nsamps_remain, replace=False)
	nsamps_per_model[rand_model_idx] = nsamps_per_model[rand_model_idx] + 1

	# Define the number of SMB samples to produce
	nmsamps = int(np.ceil(np.sqrt(nsamps)))
	ntsamps = nmsamps

	# Generate perfectly correlated samples in time for SMB
	z=rng.standard_normal(ntsamps)[:,np.newaxis]
	zit = inttemp_mean + (inttemp_sd * z)

	# Correlation between antsmb and antdyn
	#fraction=rng.random(nmsamps * ntsamps)
	fraction = None

	# Project the SMB portion of each ice sheet
	antsmb=project_antsmb(zit, my_fit, nmsamps, ntsamps, fraction=fraction)

	# Reduce the SMB years to just the target years
	year_idx = np.isin(data_years, years)
	data_years = data_years[year_idx]
	antsmb = antsmb[:,:,year_idx]

	# Convert the SMB values to mm and flatten
	antsmb *= 1000.0
	antsmb = antsmb.reshape(-1, antsmb.shape[-1])
	antsmb = antsmb[:nsamps,:]

	# Initialize data structure for the larmip samples
	eais_samps = np.full((nsamps, nyears), np.nan)
	wais_samps = np.full((nsamps, nyears), np.nan)
	pen_samps = np.full((nsamps, nyears), np.nan)
	ais_samps = np.full((nsamps, nyears), np.nan)

	# Start and end indices for populating data structures
	end_idx = np.cumsum(nsamps_per_model)
	start_idx = np.append(0, end_idx[:-1])

	# Reset the seed to allow consistent sampling of ice sheet samples with
	# other modules using the 2lm temperature data
	rng = np.random.default_rng(rngseed)

	# Loop over the models
	for i in np.arange(nmodels):

		# Get a list of sample indices that are safe to sample from
		good_sample_idx = np.flatnonzero(~np.isnan(eais_samples[i,:,0]))

		# Generate the sample indices
		sample_inds = rng.choice(good_sample_idx, size=nsamps_per_model[i], replace=True)

		# Append these samples to the data structures
		eais_samps[np.arange(start_idx[i], end_idx[i]),:] = eais_samples[i,sample_inds,:]
		wais_samps[np.arange(start_idx[i], end_idx[i]),:] = wais_samples[i,sample_inds,:]
		pen_samps[np.arange(start_idx[i], end_idx[i]),:] = pen_samples[i,sample_inds,:]
		ais_samps[np.arange(start_idx[i], end_idx[i]),:] = ais_samples[i,sample_inds,:]

	# If the user wants to extrapolate projections based on rates, do so here
	if cyear_start or cyear_end:
		for i in np.arange(nsamps):
			eais_samps[i,:] = ExtrapolateRate(eais_samps[i,:], years, cyear_start, cyear_end)
			wais_samps[i,:] = ExtrapolateRate(wais_samps[i,:], years, cyear_start, cyear_end)
			pen_samps[i,:] = ExtrapolateRate(pen_samps[i,:], years, cyear_start, cyear_end)

		# Reconstitute the ais samples from the components after extrapolating the rate
		ais_samps = eais_samps + wais_samps + pen_samps

	# Write the global projections to the netCDF files
	WriteNetCDF(pipeline_id, eais_samps, years, nsamps, "EAIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, wais_samps, years, nsamps, "WAIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, pen_samps, years, nsamps, "PEN", scenario, baseyear)
	WriteNetCDF(pipeline_id, ais_samps, years, nsamps, "AIS", scenario, baseyear)
	WriteNetCDF(pipeline_id, antsmb, years, nsamps, "SMB", scenario, baseyear)
	WriteNetCDF(pipeline_id, antsmb + ais_samps, years, nsamps, "TOT", scenario, baseyear)

	# Calculate the fractional contribution from each AIS component to the total AIS contributions
	# Note: Catch the warnings from these calls and ignore. Correcting these issues takes place in
	# the next block of code.
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		eais_frac = eais_samps / ais_samps
		wais_frac = wais_samps / ais_samps
		pen_frac = pen_samps / ais_samps

	# Correct the fractional matrices for potentially missing fractions
	eais_frac = fix_nan_icefrac(eais_frac)
	wais_frac = fix_nan_icefrac(wais_frac)
	pen_frac = fix_nan_icefrac(pen_frac)

	# Apply the fraction to the SMB
	eais_smb = antsmb * eais_frac
	wais_smb = antsmb * wais_frac
	pen_smb = antsmb * pen_frac

	# Add the component SMB values to the appropriate samples.  These get passed to the localization process.
	eais_samps += eais_smb
	wais_samps += wais_smb
	pen_samps += pen_smb
	ais_samps += antsmb

	# Store the variables in a pickle for the next stage
	output = {'eais_samps': eais_samps, 'wais_samps': wais_samps, 'pen_samps': pen_samps, \
				'years': years, 'scenario': scenario, 'baseyear': baseyear}
	outfilename = "{}_projections.pkl".format(pipeline_id)
	outfile = open(os.path.join(os.path.dirname(__file__), outfilename), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(0)



def ExtrapolateRate(sample, targyears, cyear_start, cyear_end):

	# If only one of the constant rate years is provided, imply the other
	if cyear_start and not cyear_end:
		cyear_end = cyear_start + 20
	if cyear_end and not cyear_start:
		cyear_start = cyear_end - 20

	# Find the start and end projection values for the rate calculation
	proj_start = np.interp(cyear_start, targyears, sample)
	proj_end = np.interp(cyear_end, targyears, sample)

	# Calculate the rate
	rate = (proj_end - proj_start) / (cyear_end - cyear_start)

	# Make a new projection
	ext_sample = sample
	ext_sample[targyears >= cyear_end] = proj_end + (rate * (targyears[targyears >= cyear_end] - cyear_end))

	# Return this sample
	return(ext_sample)



def fix_nan_icefrac(frac):

	# Find the samples that are causing trouble
	trouble_idx = np.flatnonzero(np.isnan(frac[:,0]))

	# If there aren't any, return the original fractional matrix
	if(len(trouble_idx) == 0):
		return(frac)

	# Replace the problematic fraction time series with the mean fractional time series
	frac[trouble_idx,:] = np.nanmean(frac, axis=0)

	# Return the fixed fractional matrix
	return(frac)



'''
project_antsmb()

Projects the surface mass balance for Antarctica using the AR5 method

'''
def project_antsmb(zit, fit_dict, nr, nt, fraction=None):
	# Return projection of Antarctic SMB contribution as a cf.Field
	# zit -- cf.Field, ensemble of time-integral temperature anomaly timeseries
	# template -- cf.Field with the required shape of the output
	# fraction -- array-like, random numbers for the SMB-dynamic feedback

	# Extract relevant parameters from the fit dictionary
	pcoK = fit_dict['pcoK']
	KoKg = fit_dict['KoKg']
	mSLEoGt = fit_dict['mSLEoGt']
	smax = fit_dict['smax']

	# Generate a distribution of products of the above two factors
	pcoKg = (pcoK[0]+rng.standard_normal([nr,nt,1])*pcoK[1])*\
		(KoKg[0]+rng.standard_normal([nr,nt,1])*KoKg[1])
	meansmb = 1923 # model-mean time-mean 1979-2010 Gt yr-1 from 13.3.3.2
	moaoKg = -pcoKg * 1e-2 * meansmb * mSLEoGt # m yr-1 of SLE per K of global warming

	if fraction is None:
		fraction=rng.random([nr,nt,1])
	elif fraction.size!=nr*nt:
		raise ProjectionError('fraction is the wrong size')
	else:
		fraction.shape=(nr,nt,1)


	ainterfactor = 1 - fraction * smax

	antsmb = moaoKg * ainterfactor * zit.reshape(1,nt,-1)

	return antsmb


'''
##########################################################################################
AR5 ice sheet projection code

Note: Leaving this here for later in case I need to break down the AIS SMB into east
and west components.  Currently, the ice fraction input file assumes total AIS (smb+dyn).
At the moment, I'm mostly concerned with generating global projections from this module,
so I'll leave this all for another day.

##########################################################################################
def ar5_project_icesheets(rng_seed, nmsamps, ntsamps, nsamps, pipeline_id):

	# Load in the icesheet fraction data-------------------------------------------
	# Note: These were derived from the Kopp14 workflow.

	# Initialize the data structures
	ice_frac = []
	ice_region_names = []

	# Open the icesheet fraction file
	ice_frac_file = os.path.join(os.path.dirname(__file__), "icesheet_fraction.txt")
	with open(ice_frac_file, 'r') as fp:

		# Get the fraction years from the header line
		header_items = re.split(",\s*", fp.readline())
		ice_frac_years = np.array([int(x) for x in header_items[1:]])

		# Read in the rest of the files
		for line in fp:
			line = line.rstrip()

			# Split the line into the region name and the fractions then append to data structures
			line_parts = re.split(",\s*", line)
			ice_region_names.append(line_parts[0])
			ice_frac.append([float(x) for x in line_parts[1:]])

	# Convert the fraction data structure into a numpy array
	ice_frac = np.array(ice_frac)

	# Subset the fraction data to the years of interest
	year_idx = np.isin(ice_frac_years, data_years)
	ice_frac = ice_frac[:,year_idx]

	# Reshape the samples and fraction data structures for broadcasting
	antnet = antnet[:,np.newaxis,:]
	ice_frac = ice_frac[np.newaxis,:,:]

	# Apply the regional fractions to the global projections
	aissamps = antnet * ice_frac

	# Save the global glacier and ice caps projections to a pickle
	output = {"gissamps": greennet, "aissamps": antnet, "totsamps": totalnet, \
		"waissamps": aissamps[:,0,:], "eaissamps": aissamps[:,1,:], "data_years": data_years}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	# Write the estimated east and west AIS contributions to netCDF files
	WriteNetCDF(aissamps[:,0,:], "WAIS", data_years, scenario, pipeline_id)
	WriteNetCDF(aissamps[:,1,:], "EAIS", data_years, scenario, pipeline_id)

	return(0)
'''



def WriteNetCDF(pipeline_id, global_samps, years, nsamps, ice_source, scenario, baseyear):

	# Write the total global projections to a netcdf file
	nc_filename = os.path.join(os.path.dirname(__file__), "{0}_{1}_globalsl.nc".format(pipeline_id, ice_source))
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(years))
	samp_dim = rootgrp.createDimension("samples", nsamps)
	loc_dim = rootgrp.createDimension("locations", 1)

	# Populate dimension variables
	year_var = rootgrp.createVariable("years", "i4", ("years",))
	samp_var = rootgrp.createVariable("samples", "i8", ("samples",))
	loc_var = rootgrp.createVariable("locations", "i8", ("locations",))
	lat_var = rootgrp.createVariable("lat", "f4", ("locations",))
	lon_var = rootgrp.createVariable("lon", "f4", ("locations",))

	# Create a data variable
	samps = rootgrp.createVariable("sea_level_change", "f4", ("samples", "years", "locations"), zlib=True, complevel=4)

	# Assign attributes
	rootgrp.description = "Global SLR contribution from {0} from the LARMIP2 workflow".format(ice_source)
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}".format(pipeline_id)
	rootgrp.scenario = scenario
	rootgrp.baseyear = baseyear
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = years
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = global_samps[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	# Close the netcdf
	rootgrp.close()

	# Done
	return(None)


if __name__ == '__main__':

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the IPCC AR6 LARMIP ice sheet projection stage.",\
	epilog="Note: This is meant to be run as part of the ipccar6 module set within the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to draw (default = 100)", default=100, type=int)
	parser.add_argument('--replace', help="Allow sampling with replacement (default = 1)", choices=(0,1), type=int, default=1)
	parser.add_argument('--seed', help="Seed for the random number generator (default = 1234)", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
	parser.add_argument('--cyear_start', help="Constant rate calculation for projections starts at this year", default=None, type=int)
	parser.add_argument('--cyear_end', help="Constant rate calculation for projections ends at this year", default=None, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Run the projection stage with the provided arguments
	ipccar6_project_larmipicesheet(args.nsamps, args.pipeline_id, args.cyear_start, args.cyear_end, args.replace, args.seed)

	exit()