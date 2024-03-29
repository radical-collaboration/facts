#!/bin/sh

# ------------------------------------------------------------------------------

pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib h5py memory-profiler
tar -xvf ./modules-data/tlm_sterodynamics_preprocess_data.tgz
tar -xvf ./modules-data/tlm_sterodynamics_cmip6_data.tgz
tar -xvf ./modules-data/ipccar6_climate_data.tgz

python3 tlm_sterodynamics_preprocess.py "--scenario" "ssp585" "--baseyear" "2005" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10" "--pipeline_id" "configTest.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics preprocess complete ==="

python3 tlm_sterodynamics_fit.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics fit complete ==="

python3 tlm_sterodynamics_project.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics" "--nsamps" "2000"

echo "=== tlm/sterodynamics project complete ==="

python3 tlm_sterodynamics_postprocess.py "--nsamps" "2000" "--pipeline_id" "capstone.modules.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics postprocess complete ==="
echo "=== tlm/sterodynamics complete ==="