import numpy as np
from netCDF4 import Dataset
import argparse
import h5py
import time
import os
import sys
import pickle
import fnmatch
import re

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


def ReadResponseFunctions(model, Tlen=None):

	# Read in the RF from the files
	fname = "./RFunctions/RF_{}_BM08_R1.dat".format(model)
	with open(fname) as f:
		r1 = np.array([float(row) for row in f])

	fname = "./RFunctions/RF_{}_BM08_R2.dat".format(model)
	with open(fname) as f:
		r2 = np.array([float(row) for row in f])

	fname = "./RFunctions/RF_{}_BM08_R3.dat".format(model)
	with open(fname) as f:
		r3 = np.array([float(row) for row in f])

	fname = "./RFunctions/RF_{}_BM08_R4.dat".format(model)
	with open(fname) as f:
		r4 = np.array([float(row) for row in f])

	fname = "./RFunctions/RF_{}_BM08_R5.dat".format(model)
	with open(fname) as f:
		r5 = np.array([float(row) for row in f])

	# Pad with zeros if necessary
	if Tlen is not None:

		zerofill = np.zeros([Tlen-200])
		r1 = np.concatenate((r1,zerofill))
		r2 = np.concatenate((r2,zerofill))
		r3 = np.concatenate((r3,zerofill))
		r4 = np.concatenate((r4,zerofill))
		r5 = np.concatenate((r5,zerofill))

	# Done
	return(r1,r2,r3,r4,r5)



def GenerateSample(M, RF, Tlen):

	# Initialize the return variable
	sample = [0]

	# Loop over the time steps
	for t in range(1, Tlen):

		# Calculate the change in sea-level
		dsl = np.sum(M[:t]*RF[t:0:-1])

		# Append this to the output array
		sample.append(dsl)

	return(np.array(sample))


def GetModelName(fname):

	model_name = None
	match = re.search(r"^RF_(\w+)_BM08_R1\.dat", fname)
	if match is not None:
		model_name = match.group(1)

	return(model_name)


def GetAvailableModels(model_dir):

	# Initialize return list
	modellist = []

	# Get a list of available R function files
	model_files = os.listdir(model_dir)

	# Generate the model list
	for this_file in model_files:
		model_name = GetModelName(this_file)
		if model_name is None:
			continue
		modellist.append(model_name)

	# Return the list
	return(modellist)



def larmip_project_icesheet(pipeline_id, nsamps, targyears, baseyear, seed, models, cyear_start, cyear_end):

	# Load the preprocessed and calibration data
	preprocess_file = "{}_preprocess.pkl".format(pipeline_id)
	with open(preprocess_file, 'rb') as f:
		preprocess_data = pickle.load(f)

	fit_file = "{}_fit.pkl".format(pipeline_id)
	with open(fit_file, 'rb') as f:
		fit_data = pickle.load(f)

	smb_file = "{}_fitsmb.pkl".format(pipeline_id)
	with open(smb_file, 'rb') as f:
		smb_data = pickle.load(f)

	# Extract the preprocessed and calibration data
	SAT = preprocess_data["SAT"]
	Time = preprocess_data["Time"]
	NumTensemble = preprocess_data["NumTensemble"]
	Tlen = preprocess_data["Tlen"]
	scenario = preprocess_data["scenario"]

	NumOmodel = fit_data["NumOmodel"]
	OS_NoDelay_R1 = fit_data["OS_NoDelay_R1"]
	OS_NoDelay_R2 = fit_data["OS_NoDelay_R2"]
	OS_NoDelay_R3 = fit_data["OS_NoDelay_R3"]
	OS_NoDelay_R4 = fit_data["OS_NoDelay_R4"]
	OS_WiDelay_R1 = fit_data["OS_WiDelay_R1"]
	OS_WiDelay_R2 = fit_data["OS_WiDelay_R2"]
	OS_WiDelay_R3 = fit_data["OS_WiDelay_R3"]
	OS_WiDelay_R4 = fit_data["OS_WiDelay_R4"]
	OS_Delay_R1 = fit_data["OS_Delay_R1"]
	OS_Delay_R2 = fit_data["OS_Delay_R2"]
	OS_Delay_R3 = fit_data["OS_Delay_R3"]
	OS_Delay_R4 = fit_data["OS_Delay_R4"]
	MeltSensitivity = fit_data["MeltSensitivity"]

	# Available models
	#available_models = ["AISM_VUB", "BISI_LBL", "CISM_NCA", "FETI_ULB", "GRIS_LSC",
	#"IMAU_UU", "ISSM_JPL", "ISSM_UCI", "MALI_DOE", "PISM_AWI", "PISM_DMI",
	#"PISM_PIK", "PISM_VUW", "PS3D_PSU", "SICO_ILTS", "UA_UNN"]
	available_models = GetAvailableModels("./RFunctions")

	if models is None:
		models = available_models

	# Number of models
	nmodels = len(models)

	# Set the rng seed
	rng = np.random.default_rng(seed)

	# How many samples per model?
	samps_per_model = np.array([nsamps // nmodels for x in range(nmodels)])
	remainder_samps = nsamps % nmodels
	samps_per_model[rng.choice(nmodels, size=remainder_samps, replace=False)] += 1
	rng.shuffle(samps_per_model)

	# Initialize the sea-level sample variables
	sl_r1 = []	# EAIS
	sl_r2 = []	# Ross
	sl_r3 = []	# Amundsen
	sl_r4 = []	# Weddell
	sl_r5 = []	# Peninsula
	sl_smb = []	# Surface Mass Balance over all Antarctica

	# Set the rng seed again to help with diagnostics
	rng = np.random.default_rng(seed)

	# Loop over the requested models
	tempcount = 0;
	for model_idx, this_model in enumerate(models):

		# Read in the appropriate model's response functions
		RF_R1, RF_R2, RF_R3, RF_R4, RF_R5 = ReadResponseFunctions(this_model, Tlen)

		# Loop over the number of samples for this model
		for i in np.arange(samps_per_model[model_idx]):

			# Choose a random forcing from the temperature data
			temp_idx = tempcount
			tempcount = tempcount + 1
			Temp = np.array(SAT[:,temp_idx])

			# Choose a random ocean model
			ocean_model_idx = rng.integers(0,NumOmodel-1)

			OS_R1 = OS_WiDelay_R1[ocean_model_idx]
			OS_R2 = OS_WiDelay_R2[ocean_model_idx]
			OS_R3 = OS_WiDelay_R3[ocean_model_idx]
			OS_R4 = OS_WiDelay_R4[ocean_model_idx]
			OS_R5 = OS_WiDelay_R4[ocean_model_idx]

			tau_R1 = int(OS_Delay_R1[ocean_model_idx])
			tau_R2 = int(OS_Delay_R2[ocean_model_idx])
			tau_R3 = int(OS_Delay_R3[ocean_model_idx])
			tau_R4 = int(OS_Delay_R4[ocean_model_idx])
			tau_R5 = int(OS_Delay_R4[ocean_model_idx])

			# Apply the delay to the temperature data
			Temp_R1 = np.append(np.zeros(tau_R1),Temp[:Tlen-tau_R1])
			Temp_R2 = np.append(np.zeros(tau_R2),Temp[:Tlen-tau_R2])
			Temp_R3 = np.append(np.zeros(tau_R3),Temp[:Tlen-tau_R3])
			Temp_R4 = np.append(np.zeros(tau_R4),Temp[:Tlen-tau_R4])
			Temp_R5 = np.append(np.zeros(tau_R5),Temp[:Tlen-tau_R5])

			# select melting sensitivity
			MS_R1 = rng.uniform(MeltSensitivity[0],MeltSensitivity[1])
			MS_R2 = rng.uniform(MeltSensitivity[0],MeltSensitivity[1])
			MS_R3 = rng.uniform(MeltSensitivity[0],MeltSensitivity[1])
			MS_R4 = rng.uniform(MeltSensitivity[0],MeltSensitivity[1])
			MS_R5 = rng.uniform(MeltSensitivity[0],MeltSensitivity[1])

			# Compose forcing time series
			M_R1 = MS_R1*OS_R1*Temp_R1
			M_R2 = MS_R2*OS_R2*Temp_R2
			M_R3 = MS_R3*OS_R3*Temp_R3
			M_R4 = MS_R4*OS_R4*Temp_R4
			M_R5 = MS_R5*OS_R5*Temp_R5

			M_R1[M_R1 < 0.0] = 0.0
			M_R2[M_R2 < 0.0] = 0.0
			M_R3[M_R3 < 0.0] = 0.0
			M_R4[M_R4 < 0.0] = 0.0
			M_R5[M_R5 < 0.0] = 0.0

			# Generate an ensemble projection
			this_sample = GenerateSample(M_R1, RF_R1, Tlen)
			sl_r1.append(this_sample)
			this_sample = GenerateSample(M_R2, RF_R2, Tlen)
			sl_r2.append(this_sample)
			this_sample = GenerateSample(M_R3, RF_R3, Tlen)
			sl_r3.append(this_sample)
			this_sample = GenerateSample(M_R4, RF_R4, Tlen)
			sl_r4.append(this_sample)
			this_sample = GenerateSample(M_R5, RF_R5, Tlen)
			sl_r5.append(this_sample)

			# Generate a sample of SMB
			this_sample = project_antsmb(Temp, Time, baseyear,rng, smb_data)
			sl_smb.append(this_sample)


	# Convert to numpy arrays and from meters to millimeters
	sl_r1 = np.array(sl_r1) * 1000.0
	sl_r2 = np.array(sl_r2) * 1000.0
	sl_r3 = np.array(sl_r3) * 1000.0
	sl_r4 = np.array(sl_r4) * 1000.0
	sl_r5 = np.array(sl_r5) * 1000.0
	sl_smb = np.array(sl_smb) * 1000.0

	# Center the projections to the baseyear
	baseyear_idx = np.flatnonzero(Time == baseyear)
	ref_vals = sl_r1[:,baseyear_idx]
	sl_r1 = sl_r1 - ref_vals
	ref_vals = sl_r2[:,baseyear_idx]
	sl_r2 = sl_r2 - ref_vals
	ref_vals = sl_r3[:,baseyear_idx]
	sl_r3 = sl_r3 - ref_vals
	ref_vals = sl_r4[:,baseyear_idx]
	sl_r4 = sl_r4 - ref_vals
	ref_vals = sl_r5[:,baseyear_idx]
	sl_r5 = sl_r5 - ref_vals
	ref_vals = sl_smb[:,baseyear_idx]
	sl_smb = sl_smb - ref_vals

	# Subset the projections for particular target years
	targyear_idx = np.array([x for x in np.arange(len(Time)) if Time[x] in targyears])
	sl_r1 = sl_r1[:,targyear_idx]
	sl_r2 = sl_r2[:,targyear_idx]
	sl_r3 = sl_r3[:,targyear_idx]
	sl_r4 = sl_r4[:,targyear_idx]
	sl_r5 = sl_r5[:,targyear_idx]
	sl_smb = sl_smb[:,targyear_idx]

	if cyear_start or cyear_end:
			for i in np.arange(nsamps):
				sl_r1[i,:] = ExtrapolateRate(sl_r1[i,:], targyears, cyear_start, cyear_end)
				sl_r2[i,:] = ExtrapolateRate(sl_r2[i,:], targyears, cyear_start, cyear_end)
				sl_r3[i,:] = ExtrapolateRate(sl_r3[i,:], targyears, cyear_start, cyear_end)
				sl_r4[i,:] = ExtrapolateRate(sl_r4[i,:], targyears, cyear_start, cyear_end)
				sl_r5[i,:] = ExtrapolateRate(sl_r5[i,:], targyears, cyear_start, cyear_end)
				sl_smb[i,:] = ExtrapolateRate(sl_smb[i,:], targyears, cyear_start, cyear_end)

	# Combine all Antarctic regions into Antarctica projection
	sl_su = sl_r1 + sl_r2 + sl_r3 + sl_r4 + sl_r5

	# Sum up the WAIS components
	wais_samples = sl_r2 + sl_r3 + sl_r4

	# Write the global projections to netcdf files
	WriteNetCDF(sl_su + sl_smb, None, targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(sl_r1, "EAIS", targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(sl_r5, "PEN", targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(wais_samples, "WAIS", targyears, baseyear, scenario, nsamps, pipeline_id)
	WriteNetCDF(sl_smb, "SMB", targyears, baseyear, scenario, nsamps, pipeline_id)

	# Distribute the SMB across the AIS regions
	r1_frac = sl_r1 / sl_su
	r2_frac = sl_r2 / sl_su
	r3_frac = sl_r3 / sl_su
	r4_frac = sl_r4 / sl_su
	r5_frac = sl_r5 / sl_su

	r1_smb = r1_frac * sl_smb
	r2_smb = r2_frac * sl_smb
	r3_smb = r3_frac * sl_smb
	r4_smb = r4_frac * sl_smb
	r5_smb = r5_frac * sl_smb

	sl_r1 += r1_smb
	sl_r2 += r2_smb
	sl_r3 += r3_smb
	sl_r4 += r4_smb
	sl_r4 += r5_smb

	# Save the projections to a pickle
	output = {"sl_r1": sl_r1, "sl_r2": sl_r2, "sl_r3": sl_r3, "sl_r4": sl_r4, "sl_r5": sl_r5, \
				"targyears": targyears, "baseyear": baseyear, "models": models, "scenario": scenario}
	outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
	pickle.dump(output, outfile)
	outfile.close()

	return(None)


'''
project_antsmb()

Projects the surface mass balance for Antarctica using the AR5 method.
Adapted from code provided by Jonathan Gregory.

'''
def project_antsmb(temp_sample, years, baseyear, rng,fit_dict,):

	# Initialize the smb values for this sample
	antsmb = np.zeros(len(temp_sample))

	# Begin integrating temperature from the base year
	baseyear_idx = int(np.flatnonzero(years == baseyear))
	int_temp = np.cumsum(temp_sample[baseyear_idx:])

	# Extract relevant parameters from the fit dictionary
	pcoK = fit_dict['pcoK']
	KoKg = fit_dict['KoKg']
	mSLEoGt = fit_dict['mSLEoGt']
	smax = fit_dict['smax']

	# Generate a distribution of products of the above two factors
	#pcoKg = (pcoK[0]+rng.standard_normal([nr,nt,1])*pcoK[1])*\
	#	(KoKg[0]+rng.standard_normal([nr,nt,1])*KoKg[1])
	pcoKg = (pcoK[0]+rng.standard_normal(1)*pcoK[1])*\
		(KoKg[0]+rng.standard_normal(1)*KoKg[1])
	meansmb = 1923 # model-mean time-mean 1979-2010 Gt yr-1 from 13.3.3.2
	moaoKg = -pcoKg * 1e-2 * meansmb * mSLEoGt # m yr-1 of SLE per K of global warming

	antsmb[baseyear_idx:] = moaoKg * (1-smax) * int_temp

	return antsmb



def WriteNetCDF(slr, region, targyears, baseyear, scenario, nsamps, pipeline_id):

	# Write the total global projections to a netcdf file
	if region is None:
		nc_filename = os.path.join(os.path.dirname(__file__), "{}_globalsl.nc".format(pipeline_id))
		nc_description = "Global SLR contribution from Antarctica using the LARMIP module"
	else:
		nc_filename = os.path.join(os.path.dirname(__file__), "{}_{}_globalsl.nc".format(pipeline_id, region))
		nc_description = "Global SLR contribution from Antarctica ({}) using the LARMIP module".format(region)
	rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

	# Define Dimensions
	year_dim = rootgrp.createDimension("years", len(targyears))
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
	rootgrp.description = nc_description
	rootgrp.history = "Created " + time.ctime(time.time())
	rootgrp.source = "FACTS: {0}. ".format(pipeline_id)
	rootgrp.baseyear = baseyear
	rootgrp.scenario = scenario
	samps.units = "mm"

	# Put the data into the netcdf variables
	year_var[:] = targyears
	samp_var[:] = np.arange(nsamps)
	samps[:,:,:] = slr[:,:,np.newaxis]
	lat_var[:] = np.inf
	lon_var[:] = np.inf
	loc_var[:] = -1

	return(None)





if __name__ == "__main__":

	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the projection stage for the LARMIP ice sheet module",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

	# Define the command line arguments to be expected
	parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module", required=True)
	parser.add_argument('--models', help="List of model names to include in the sampling process", nargs="+", default=None)
	parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
	parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)
	parser.add_argument('--cyear_start', help="Constant rate calculation for projections starts at this year", default=None, type=int)
	parser.add_argument('--cyear_end', help="Constant rate calculation for projections ends at this year", default=None, type=int)

	# Parse the arguments
	args = parser.parse_args()

	# Run the larmip code
	targyears = np.arange(args.pyear_start, args.pyear_end+1, args.pyear_step)
	larmip_project_icesheet(nsamps=args.nsamps, pipeline_id=args.pipeline_id, cyear_start=args.cyear_start, cyear_end=args.cyear_end,targyears=targyears, seed=args.seed, baseyear=args.baseyear, models=args.models)


	sys.exit()