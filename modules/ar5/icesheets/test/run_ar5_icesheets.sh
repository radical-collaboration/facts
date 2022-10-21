#!/bin/bash
SCRIPT_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}") 
cd $SCRIPT_RELATIVE_DIR
SCRIPT_DIR=`pwd`

WORKDIR=/scratch/`whoami`/test.`date +%s`
mkdir $WORKDIR
cd $WORKDIR

pwd

cp -L -r $SCRIPT_DIR/../* .
LFILE=test/location_list.lst

for i in data/*
do
    tar xzvf $i  2>&1 | grep -v 'Ignoring'
done



python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585 --scenario ssp585 --baseyear 2005 --tlm_data 1
python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585
python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --chunksize 662 --locationfile ${LFILE}

mkdir $SCRIPT_DIR/output_global
mkdir $SCRIPT_DIR/output_local
mv *globalsl* $SCRIPT_DIR/output_global
mv *localsl* $SCRIPT_DIR/output_local

