# GMSLR projection program used for IPCC WG1 AR5
# Translated from IDL to Python 2.7 by Jonathan Gregory 23.10.19

import os
import numpy as np
import netCDF4 as nc
import pandas as pd
import time
from netCDF4 import Dataset
import pickle
import sys
import re
import argparse

class ProjectionError(Exception):
	pass
	

def ebm3_thermalexpansion_project(scenario, climate_data_file, coef_file, params_file, pyear_start, pyear_end, pyear_step, nsamps, pipeline_id, seed, baseyear):
	# constants 
    a = 6.37*1e6
    earth_area = 4*np.pi*a**2
	
    path = os.path.dirname(__file__)

    # heat capacity of each layer from fair2. Needed to compute OHC
    fparam = pd.read_csv(params_file)
    c1 = fparam['clim_c1']
    c2 = fparam['clim_c2']
    c3 = fparam['clim_c3']
	
    nsims=1001

    rng = np.random.default_rng(seed)
    if nsamps > nsims:
        run_idx = np.arange(nsims)
        sample_idx = rng.choice(nsims, nsamps, nsamps>nsims)
    else:
        run_idx = rng.choice(nsims, nsamps, nsamps>nsims)
        sample_idx = np.arange(nsamps)
    
    c1 = c1[sample_idx]
    c2 = c2[sample_idx]
    c3 = c3[sample_idx]

    # create target years array
    targyears = np.arange(pyear_start, pyear_end+1, pyear_step)

    #temperature output from fair2. Needed t compute OHC
    ds = Dataset(climate_data_file)
    gsat = ds[scenario]['surface_temperature'][:]
    deepoceant = ds[scenario]['deep_ocean_temperature'][:]
    years = ds[scenario]['years'][:]

    #   Expansion coefficients. Needed to to convert OCH to global thermal expansion (GTE)
    ds = Dataset(coef_file)
    include_models = ds['model'][:]
    eeh3 = ds['expcoefs'][:]

    # estimate OHC
    up = gsat*np.array(c1)
    mid = deepoceant[:,:,0]*np.array(c2)
    low = deepoceant[:,:,1]*np.array(c3)
    n = np.zeros((len(years), nsamps))
    n[1:,:] = up[1:,:]-up[:-1,:] + mid[1:,:]-mid[:-1,:] + low[1:,:]-low[:-1,:]
    ohc_samps = np.cumsum(n*earth_area, 0) * 365*24*3.6*1e3*1e-24

    # Generate samples assuming normal distribution
    rng = np.random.default_rng(seed)
    expcoef_samps = rng.normal(loc=np.mean(eeh3), scale=np.std(eeh3), size=(nsamps,1))
 
    # Produce the projection samples
    gte_samps = ohc_samps * expcoef_samps.flatten()

    # Center these samples on the baseyear
    baseyear_idx = np.flatnonzero(years == baseyear)
    gte_samps = gte_samps - gte_samps[baseyear_idx,:]

    # Subset the samples for the projection years
    targyear_idx = np.isin(years, targyears)
    gte_samps = gte_samps[targyear_idx,:]

    # Invert the dimensions of the variable and convert from m to mm
    gte_samps *= 1000.
    gte_samps = gte_samps.T

    # Save the projections to a pickle
    output = {"thermsamps": gte_samps, "targyears": targyears, "baseyear": baseyear, \
 			"include_models": include_models, "scenario": scenario}
    outfile = open(os.path.join(os.path.dirname(__file__), "{}_projections.pkl".format(pipeline_id)), 'wb')
    pickle.dump(output, outfile)
    outfile.close()

    # Write the total global projections to a netcdf file
    nc_filename = os.path.join(os.path.dirname(__file__), "{0}_globalsl.nc".format(pipeline_id))
    rootgrp = Dataset(nc_filename, "w", format="NETCDF4")

    # Define Dimensions
    nyr = len(targyears)
    year_dim = rootgrp.createDimension("years", nyr)
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
    rootgrp.description = "Global SLR contribution from Thermal Expansion according to Three-Layer Model workflow"
    rootgrp.history = "Created " + time.ctime(time.time())
    rootgrp.source = "FACTS: {0}".format(pipeline_id)
    rootgrp.scenario = scenario
    rootgrp.baseyear = baseyear
    rootgrp.comment = "Included Models: " + ",".join([str(x) for x in include_models])
    samps.units = "mm"

    # Put the data into the netcdf variables
    year_var[:] = targyears
    samp_var[:] = np.arange(nsamps)
    samps[:,:,:] = gte_samps[:,:,np.newaxis]
    lat_var[:] = np.inf
    lon_var[:] = np.inf
    loc_var[:] = -1

    # Close the netcdf
    rootgrp.close()
    
    return(0)

	

if __name__ == '__main__':
	
	# Initialize the command-line argument parser
	parser = argparse.ArgumentParser(description="Run the glacier projection stage for the AR5 SLR projection workflow",\
	epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")
	
	# Define the command line arguments to be expected
	parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)", default='ssp585')
	parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
	parser.add_argument('--params_file', help='Full path to calibrated constraints params file', default='calibrated_constrained_parameters.csv')
	parser.add_argument('--coef_file', help='Full path to expansion coefficient file', default='scmpy3LM_RCMIP_CMIP6calpm_n18_expcoefs.nc')
	parser.add_argument('--nsamps', help="Number of samples to generate [default=1000]", default=1000, type=int)
	parser.add_argument('--pyear_start', help="Projection year start [default=2020]", default=2020, type=int)
	parser.add_argument('--pyear_end', help="Projection year end [default=2100]", default=2150, type=int)
	parser.add_argument('--pyear_step', help="Projection year step [default=10]", default=10, type=int)
	parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)
	parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
	parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

	
	# Parse the arguments
	args = parser.parse_args()
	
	# Run the projection process on the files specified from the command line argument
	ebm3_thermalexpansion_project(args.scenario, args.climate_data_file, args.coef_file, args.params_file, args.pyear_start, args.pyear_end, args.pyear_step, args.nsamps, args.pipeline_id, args.seed, args.baseyear)
	
	exit()
