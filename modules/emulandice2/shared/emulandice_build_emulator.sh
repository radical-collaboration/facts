#!/bin/bash

# check for results directory
#if [[ ! -d results ]]; then
#    mkdir results
#fi

# set up environment
source emulandice_environment.sh

# run emulandice
Rscript -e "library(emulandice2)" -e "source('emulandice2/emulator_build.R')" AIS 0 2300
