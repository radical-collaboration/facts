#!/bin/bash

# check for results directory
if [[ ! -d results ]]; then
    mkdir results
fi

# set up environment
source emulandice_environment.sh --Rscript

# run emulandice
emulandice_dataset=$1
nsamps=$2
ice_sources=$3
if [ "$ice_sources" == "" ]; then
    echo Rscript -e "source('packrat/init.R')" -e "library(emulandice)" -e "emulandice::main('decades',dataset='$emulandice_dataset',N_FACTS=$nsamps)"
    Rscript -e "source('packrat/init.R')" -e "library(emulandice)" -e "emulandice::main('decades',dataset='$emulandice_dataset',N_FACTS=$nsamps)"
else
    echo Rscript -e "source('packrat/init.R')" -e "library(emulandice)" -e "emulandice::main('decades',dataset='$emulandice_dataset',N_FACTS=$nsamps,ice_sources=c('$ice_sources'))"
    Rscript -e "source('packrat/init.R')" -e "library(emulandice)" -e "emulandice::main('decades',dataset='$emulandice_dataset',N_FACTS=$nsamps,ice_sources=c('$ice_sources'))"
fi

