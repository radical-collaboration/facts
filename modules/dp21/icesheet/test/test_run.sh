#!/bin/bash

python dp21_preprocess_icesheet.py --pipeline_id icesheets-dp21-icesheet-ssp585 --scenario rcp85 --baseyear 2005
python dp21_fit_icesheet.py --pipeline_id icesheets-dp21-icesheet-ssp585
python dp21_project_icesheet.py --pipeline_id icesheets-dp21-icesheet-ssp585 --nsamps 20000 --seed 4321 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
python dp21_postprocess_icesheet.py --pipeline_id icesheets-dp21-icesheet-ssp585 --chunksize 662 --locationfile $LFILE

