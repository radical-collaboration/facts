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

def ebm3_thermalexpansion_preprocess(scenario, pipeline_id, fname):

    # Constants
    a = 6.37*1e6
    earth_area = 4*np.pi*a**2
    
    # Working directory
    path = os.path.dirname(__file__)

    # Load Ocean Heat Content
    # heat capacity of each layer from fair2. Needed to compute OHC
    path2para = 'calibrated_constrained_parameters.csv'
    fparam = pd.read_csv(path2para)

    c1 = fparam['clim_c1']
    c2 = fparam['clim_c2']
    c3 = fparam['clim_c3']

    #temperature output from fair2. Needed t compute OHC
    #path = '/Users/vmalagonsantos/Library/CloudStorage/OneDrive-NIOZ/GitHub/project_tsl/fair2_output/' + scenario + '.temperature.fair.temperature_climate.nc'
    
    ds = nc.Dataset(fname)
    gsat = ds[scenario]['surface_temperature'][:]
    deepoceant = ds[scenario]['deep_ocean_temperature'][:]
    years = ds[scenario]['years'][:]

    # Expansion coefficients. Needed to to convert OCH to global thermal expansion (GTE)
    path2exp = 'scmpy3LM_RCMIP_CMIP6calpm_n18_expcoefs.nc'
    ds = nc.Dataset(path2exp)
    include_models = ds['model'][:]
    eeh3 = ds['expcoefs'][:]

    # Store preprocessed data in pickles
    output = {'ohc_samps': ohc_samps, 'expcoefs': expcoefs, \
              'expcoefs_models': expcoefs_models, \
              'data_years': data_years, 'scenario': scenario}

    # Write the configuration to a file
    outdir = os.path.dirname(__file__)
    outfile = open(os.path.join(outdir, "{}_tlmdata.pkl".format(pipeline_id)), 'wb')
    pickle.dump(output, outfile)
    outfile.close()

    return (None)

if __name__ == '__main__':

    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the pre-processing stage for the TLM ocean dynamics workflow", \
                                     epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

    # Define the command line arguments to be expected
    parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)", default='ssp585')
    parser.add_argument('--climate_data_file', help="NetCDF4/HDF5 file containing surface temperature data", type=str)
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
    
    # Parse the arguments
    args = parser.parse_args()

    # Define the names of the intermediate files this run will generate
    outdir = os.path.dirname(__file__)
    tlmfile = os.path.join(outdir, "{}_tlmdata.pkl".format(args.pipeline_id))

    # Run the 2-layer model preprocessing if the intermediate file is not present
    ebm3_thermalexpansion_preprocess(args.scenario, args.pipeline_id, args.climate_data_file)

    # Done
    sys.exit()