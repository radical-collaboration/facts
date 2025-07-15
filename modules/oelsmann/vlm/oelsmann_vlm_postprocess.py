#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:32:46 2024

Projects VLM based on Oelsmann's approach, combining:
    
    - Caron et al. (2018) statistics.
    - InSAR-bases products.
    
@author: vmalagonsantos

"""

#%% 0. Packages and functions

import xarray as xr
import numpy as np
import time
import os
import sys
import argparse

def oelmanns_vlm_postprocess(pipeline_id, nsamps, seed, pyear_start, pyear_end, pyear_step, locationfile, baseyear, vlmdir):

    projyears = np.arange(pyear_start, pyear_end+1, pyear_step)
    filenames = [vlmdir + f for f in os.listdir(vlmdir) if not f.startswith('.')]
    
    nsims = len(filenames)
    rng = np.random.default_rng(seed)
    if nsamps > nsims:
        run_idx = np.arange(nsims)
        sample_idx = rng.choice(nsims, nsamps, nsamps>nsims)
    else:
        run_idx = rng.choice(nsims, nsamps, nsamps>nsims)
        sample_idx = np.arange(nsamps)
        print(sample_idx)
        
    filenames = [filenames[f] for f in sample_idx]  
    
    # ds = xr.open_mfdataset(filenames)
    ds = xr.open_mfdataset(filenames, combine="nested", concat_dim="samples")
    
    # invert values
    ds['sea_level_change'] = ds['sea_level_change'] * -1
    
    vlm_out = ds[["sea_level_change"]]
    # Add "lat" and "lon" as data variable in output, pulling values from the first file.
    vlm_out["lat"] = ds["lat"].isel(samples=0).values
    vlm_out["lon"] = ds["lon"].isel(samples=0).values
    
    #% select proj years and apply baseyear
    
    vlm_out = vlm_out.sel(years=projyears) - vlm_out.sel(years=baseyear)
    vlm_out = vlm_out.chunk(dict(samples=-1))
    
    lat = ds["lat"].isel(samples=0).values
    lon = ds["lon"].isel(samples=0).values
    vlm = vlm_out['sea_level_change'].values
    years = vlm_out['years'].values
    locations = vlm_out['locations'].values
    nc_missing_value = np.nan

    nc_missing_value = np.nan
    # Generate the output xarray

    ncvar_attributes = {"description": "Local SLR contributions from VLM. Computed by Oelsmann",
            "history": "Created " + time.ctime(time.time()),
            "source": "SLR Framework: VLM Oelmanns",
            "baseyear": baseyear}

    local_out = xr.Dataset({"sea_level_change": (("samples", "years", "locations"), vlm, {"units":"mm", "missing_value":nc_missing_value}),
                            "lat": (("locations"), lat),
                            "lon": (("locations"), lon)},
                            coords={"years": years, "locations": locations, "samples": np.arange(nsamps)}, attrs=ncvar_attributes)
        # Write these samples to a temporary netcdf file
    local_out.to_netcdf("{0}_localsl.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
    local_outq = local_out.quantile([0.01,0.05,0.17,0.50,0.83,0.95,0.99], dim='samples')
    local_outq.to_netcdf("{0}_quantiles.nc".format(pipeline_id), encoding={"sea_level_change": {"dtype": "f4", "zlib": True, "complevel":4, "_FillValue": nc_missing_value}})
        

if __name__ == '__main__':

    # Initialize the command-line argument parser
    parser = argparse.ArgumentParser(description="Run the post-processing stage for Oelmanss VLM workflow",\
    epilog="Note: This is meant to be run as part of the Framework for the Assessment of Changes To Sea-level (FACTS)")

    # Define the command line arguments to be expected
    parser.add_argument('--nsamps', help="Number of samples to generate [default=20000]", default=20000, type=int)
    parser.add_argument('--seed', help="Seed value for random number generator [default=1234]", default=1234, type=int)
    parser.add_argument('--pyear_start', help="Year for which projections start [default=2000]", default=2000, type=int)
    parser.add_argument('--pyear_end', help="Year for which projections end [default=2300]", default=2300, type=int)
    parser.add_argument('--pyear_step', help="Step size in years between pyear_start and pyear_end at which projections are produced [default=5]", default=5, type=int)
    parser.add_argument('--locationfile', help="File that contains name, id, lat, and lon of points for localization", default="location.lst")
    parser.add_argument('--baseyear', help="Base year to which slr projections are centered", type=int, default=2005)
    parser.add_argument('--pipeline_id', help="Unique identifier for this instance of the module")
    parser.add_argument('--vlmdir',help='Path to VLM samples directory', default='vlm_samples/')

    # Parse the arguments
    args = parser.parse_args()

    oelmanns_vlm_postprocess(args.pipeline_id, 
                                      args.nsamps, 
                                      args.seed, 
                                      args.pyear_start, 
                                      args.pyear_end, 
                                      args.pyear_step, 
                                      args.locationfile, 
                                      args.baseyear, 
                                      args.vlmdir)

    # Done
    sys.exit()
