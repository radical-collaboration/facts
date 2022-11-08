#!/bin/bash
if [[ ! -d results ]]; then
    mkdir results
fi

# run emulandice
emulandice_dataset=$1
nsamps=$2
ice_sources=$3
Rscript -e "source('packrat/init.R')" -e "library(emulandice)" -e "emulandice::main('decades',dataset='$emulandice_dataset',N_FACTS=$nsamps)"
