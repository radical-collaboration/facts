#!/bin/bash

python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --scenario ssp585 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585
python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --nsamps 20000 --seed 4321
python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --chunksize 662 --locationfile $LFILE

