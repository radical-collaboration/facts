#!/bin/bash

# check for results directory
#if [[ ! -d results ]]; then
#    mkdir results
#fi

# set up environment
source emulandice_environment.sh

# run emulandice

ice_source=$1 # Ice source: GIS, AIS or GLA
region=$2 # Ice source region: ALL for GIS/AIS and RGI01-RGI19 for GLA
emu_file=$3 # models emulator file
climate_data_file=$4 # e.g. emulandice.ssp585.temperature.fair.temperature_climate.nc
scenario=$5 # e.g. ssp585 [could extract from filename instead?]

Rscript -e "library(emulandice2)" -e "source('emulandice_steer.R')" $ice_source $region $emu_name $climate_data_file $scenario
