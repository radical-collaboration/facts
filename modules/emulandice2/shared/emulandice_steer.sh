#!/bin/bash

# check for results directory
#if [[ ! -d results ]]; then
#    mkdir results
#fi

# set up environment
source emulandice_environment.sh


ice_source=$1 # Ice source: GIS, AIS or GLA
region=$2 # Ice source region: ALL for GIS/AIS and RGI01-RGI19 for GLA
emu_file=$3 # full path to build file: must match two settings above
climate_data_file=$4 # e.g. emulandice.ssp585.temperature.fair.temperature_climate.nc
scenario=$5 # e.g. ssp585 [could extract from filename instead?]
outdir=$6 # name of directory to put outputs - will attempt to create if it doesn't exist
seed=$7
pipeline_id=$8

if [ $# != 8 ] && [ $# != 0 ];
then
      echo "Wrong number of arguments given"
      echo "Required: ice_source region emu_file climate_data_file scenario outdir seed pipeline_id"
      echo "Or if no arguments are set, will run in test mode"
     exit 1
fi

Rscript -e "source('.Rprofile')" -e "library(emulandice2)" -e "source('emulandice_steer.R')" $ice_source $region $emu_file $climate_data_file $scenario $outdir $seed $pipeline_id

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
