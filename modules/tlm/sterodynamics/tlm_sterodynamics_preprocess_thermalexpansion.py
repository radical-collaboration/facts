import numpy as np
import pickle
import os
import sys
import re
import argparse
import fnmatch
# from DriftCorr import DriftCorr

from Import2lmData import *

''' tlm_preprocess_oceandynamics.py

This runs the preprocessing stage for the ocean dynamics component of the IPCC AR6
workflow.

Parameters:
ssp_scenario = SSP scenario of interest
modeldir = Directory that contains the ZOS and ZOSTOGA CMIP6 model data
driftcorr = Apply the drift correction?
locationfilename = File that contains points for localization
pipeline_id = Unique identifier for the pipeline running this code


'''
def tlm_preprocess_thermalexpansion(scenario, pipeline_id, fname):
    # Working directory
    path = os.path.dirname(__file__)

    # Load the ocean heat content
    ohc_dict = Import2lmData("ocean_heat_content", scenario, path, climate_fname=fname)

    # Extract the ohc samples
    ohc_samps = ohc_dict["samples"] * 1e-24
    data_years = ohc_dict["years"]

    # Load expansion coefficients
    nc = Dataset(os.path.join(path, 'scmpy2LM_RCMIP_CMIP6calpm_n18_expcoefs.nc'), 'r')
    expcoefs = nc.variables['expcoefs'][...]
    expcoefs_models = nc.variables['model'][...]
    nc.close()

    # Load GSAT RMSEs
    nc = Dataset(os.path.join(path, 'scmpy2LM_RCMIP_CMIP6calpm_n17_gsat_rmse.nc'), 'r')
    rmses = nc.variables['gsat_rmse'][...]
    rmses_models = nc.variables['model'][...]
    nc.close()

    # Store preprocessed data in pickles
    output = {'ohc_samps': ohc_samps, 'expcoefs': expcoefs, 'rmses': rmses, \
              'expcoefs_models': expcoefs_models, 'rmses_models': rmses_models, \
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
    parser.add_argument('--scenario', help="SSP scenario (i.e ssp585) or temperature target (i.e. tlim2.0win0.25)",
                        default='ssp585')

    parser.add_argument('--no_drift_corr', help="Do not apply the drift correction", type=int, choices=[0, 1],
                        default=0)

    parser.add_argument('--no_correlation', help="Do not apply the correlation between ZOS and ZOSTOGA fields",
                        type=int, choices=[0, 1], default=0)

    parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2020, type=int)
    parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)

    parser.add_argument('--pyear_step',
                        help="Step size in years between pyear_start and pyear_end at which projections are produced [default=10]",
                        default=10, type=int)

    parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization",
                        default="location.lst")

    parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2000)

    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")

    # Parse the arguments
    args = parser.parse_args()

    # Define the names of the intermediate files this run will generate
    outdir = os.path.dirname(__file__)
    configfile = os.path.join(outdir, "{}_config.pkl".format(args.pipeline_id))
    zostogafile = os.path.join(outdir, "{}_ZOSTOGA.pkl".format(args.pipeline_id))
    zosfile = os.path.join(outdir, "{}_ZOS.pkl".format(args.pipeline_id))
    tlmfile = os.path.join(outdir, "{}_tlmdata.pkl".format(args.pipeline_id))


    # Run the 2-layer model preprocessing if the intermediate file is not present
    if os.path.isfile(tlmfile):
        print("{} found, skipping TE preprocessing".format(tlmfile))
    else:
        tlm_preprocess_thermalexpansion(args.scenario, args.pipeline_id, args.climate_data_file)

    # Done
    sys.exit()
