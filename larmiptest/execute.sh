#!/bin/sh

# ------------------------------------------------------------------------------

pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib h5py memory-profiler
tar -xvf ./modules-data/ipccar6_climate_data.tgz



python3 larmip_icesheet_preprocess.py "--scenario" "ssp585" "--pipeline_id" "configTest.larmip.larmip.AIS"

echo "=== larmip preprocess step COMPLETE ==="

tar -xvf ./modules-data/larmip_icesheet_fit_data.tgz
python3 larmip_icesheet_fit.py "--pipeline_id" "configTest.larmip.larmip.AIS"

echo "=== larmip fit step COMPLETE ==="

tar -xvf ./modules-data/larmip_icesheet_project_data.tgz
python3 larmip_icesheet_project.py "--nsamps" "2000" "--baseyear" "2005" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10" "--pipeline_id" "configTest.larmip.larmip.AIS"

echo "=== larmip project step COMPLETE ==="

tar -xvf ./modules-data/grd_fingerprints_data.tgz
python3 larmip_icesheet_postprocess.py "--pipeline_id" "configTest.larmip.larmip.AIS"

echo "=== larmip postprocess step COMPLETE ==="

echo "=== larmip COMPLETE ==="