#!/bin/bash

if [ ! -d $HOME/singularity ]; then
    mkdir $HOME/singularity
fi

if [ ! -d $HOME/singularity/mongo ]; then
    mkdir $HOME/singularity/mongo
fi

if [ ! -f $HOME/singularity/mongo/singularity ]; then
    singularity build --sandbox $HOME/singularity/mongo docker://mongo
fi
