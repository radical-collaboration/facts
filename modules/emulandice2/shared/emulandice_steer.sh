#!/bin/bash

# check for results directory
#if [[ ! -d results ]]; then
#    mkdir results
#fi

# set up environment
source emulandice_environment.sh

# Run emulandice2 module

ice_source=$1 # Ice source: GIS, AIS or GLA
region=$2 # Ice source region: ALL for GIS/AIS and RGI01-RGI19 for GLA
emu_name=$3 # models_emulator_settings: e.g. "CISM_pow_exp_20", "CISM_IMAUICE_GISM_pow_exp_20"
emu_file=$4 # full path to build file: must match three settings above
climate_data_file=$5 # e.g. emulandice.ssp585.temperature.fair.temperature_climate.nc
scenario=$6 # e.g. ssp585 [could extract from filename instead?]
outdir=$7 # name of directory to put outputs - will attempt to create if it doesn't exist
seed=$8
pipeline_id=$9

echo $3
echo $4
if [ $# != 9 ] && [ $# != 0 ];
then
      echo "Insufficient arguments given"
      echo "Required: ice_source region emu_name emu_file climate_data_file scenario outdir seed pipeline_id"
      echo "If no arguments are set, will run in test mode"
     exit 1
fi

Rscript -e "library(emulandice2)" -e "source('emulandice_steer.R')" $ice_source $region $emu_name $emu_file $climate_data_file $scenario $outdir $seed $pipeline_id

# Arguments to add ____________
# Read from emulator RData file:
# baseyear # 2000
# pyear_start # 2005
# pyear_end # 2100 or 2300
# pyear_step # 5

# Other:
# nsamps # could add if want to deviate from number of FaIR samples or hard check on this

# Build options for emu_name used in PROTECT Fall 2023 meeting results:
# AIS ALL Kori_PISM_pow_exp_10
# GIS ALL CISM_IMAUICE_pow_exp_20
# GIS ALL CISM_pow_exp_20
# GLA RGI03 GloGEM_OGGM_pow_exp_20

Rscript -e "library(emulandice2)" -e "source('emulandice_steer_working.R')" $ice_source $region $emu_file $climate_data_file $scenario
