#!/bin/bash

# ------------------------------------------------------------------------------

# pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib h5py memory-profiler
tar -xvf ./modules-data/tlm_sterodynamics_preprocess_data.tgz
tar -xvf ./modules-data/tlm_sterodynamics_cmip6_data.tgz
tar -xvf ./modules-data/ipccar6_climate_data.tgz

# python3 tlm_sterodynamics_preprocess.py "--scenario" "ssp585" "--baseyear" "2005" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10" "--pipeline_id" "configTest.ocean.tlm.sterodynamics"
mprof run -o mem_tlm_preprocess.dat tlm_sterodynamics_preprocess.py "--scenario" "ssp585" "--baseyear" "2005" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10" "--pipeline_id" "configTest.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics preprocess complete ==="

# python3 tlm_sterodynamics_fit.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics"
mprof run -o mem_tlm_fit.dat tlm_sterodynamics_fit.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics fit complete ==="

# python3 tlm_sterodynamics_project.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics" "--nsamps" "2000"
mprof run -o mem_tlm_project.dat tlm_sterodynamics_project.py "--pipeline_id" "configTest.ocean.tlm.sterodynamics" "--nsamps" "2237"

echo "=== tlm/sterodynamics project complete ==="

# python3 tlm_sterodynamics_postprocess.py "--nsamps" "2000" "--pipeline_id" "capstone.modules.ocean.tlm.sterodynamics"
mprof run -o mem_tlm_postprocess.dat tlm_sterodynamics_postprocess.py "--nsamps" "2237" "--pipeline_id" "configTest.ocean.tlm.sterodynamics"

echo "=== tlm/sterodynamics postprocess complete ==="
echo "=== tlm/sterodynamics complete ==="

find . -name "*globalsl.nc" -type f -exec mv -t /opt/total_data/global/ {} +
find . -name "*localsl.nc" -type f -exec mv -t /opt/total_data/local/ {} +
# find . -name "*.nc" -type f -exec mv -t /opt/climate {} +
mv *.nc /opt/climate
