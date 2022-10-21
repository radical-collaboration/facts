#!/bin/bash

python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585 --scenario ssp585 --baseyear 2005 --tlm_data 1
python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585
python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --chunksize 662 --locationfile $LFILE
