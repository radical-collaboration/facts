#!/bin/sh

# ------------------------------------------------------------------------------

# mprof run -o mem_preprocess.dat "fair_temperature_preprocess.py" "--scenario" "ssp585" "--pipeline_id" "configTest.temperature.fair.temperature" &
tar -xvf ./modules-data/fair_temperature_preprocess_data.tgz 
        # rm /modules-data/fair_temperature_preprocess_data.tgz || rp_error pre_exec
# pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib pandas xarray dask fair==1.6.4 memory-profiler || rp_error pre_exec

echo "we lit"

python3 fair_temperature_preprocess.py "--scenario" "ssp585" "--pipeline_id" "configTest.temperature.fair.temperature"

tar -xvf ./modules-data/fair_temperature_fit_data.tgz 
python3 fair_temperature_fit.py "--pipeline_id" "configTest.temperature.fair.temperature" 


python3 fair_temperature_project.py "--pipeline_id" "configTest.temperature.fair.temperature" "--nsamps" "2000" 


python3 fair_temperature_postprocess.py "--pipeline_id" "configTest.temperature.fair.temperature" 