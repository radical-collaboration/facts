#!/bin/bash

python kopp14_preprocess_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
python kopp14_fit_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
python kopp14_project_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
python kopp14_postprocess_verticallandmotion_dev.py --pipeline_id verticallandmotion-kopp14-verticallandmotion --nsamps 20000 --seed 7331 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --chunksize 662 --locationfile $LFILE
