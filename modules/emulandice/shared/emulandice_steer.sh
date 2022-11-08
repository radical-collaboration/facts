#!/bin/bash
if [[ ! -d results ]]; then
    mkdir results
fi

# run emulandice
emulandice_dataset=$1
Rscript -e "emulandice::main('decades',dataset=\'$emulandice_dataset\')"
