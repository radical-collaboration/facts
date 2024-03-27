#!/bin/sh

# ------------------------------------------------------------------------------

pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib h5py xarray dask[array] memory-profiler

tar -xvf ./modules-data/ipccar6_climate_data.tgz
python3 FittedISMIP_GrIS_preprocess.py "--scenario" "ssp585" "--pipeline_id" "configTest.GrIS1f.FittedISMIP.GrIS"

tar -xvf ./modules-data/FittedISMIP_icesheet_fit_data.tgz
python3 FittedISMIP_GrIS_fit.py "--pipeline_id" "configTest.GrIS1f.FittedISMIP.GrIS"

python3 FittedISMIP_GrIS_project.py "--nsamps" "2000" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10" "--baseyear" "2005" "--pipeline_id" "configTest.GrIS1f.FittedISMIP.GrIS"

tar -xvf ./modules-data/grd_fingerprints_data.tgz
python3 FittedISMIP_GrIS_postprocess.py "--pipeline_id" "configTest.GrIS1f.FittedISMIP.GrIS"

ls -la
