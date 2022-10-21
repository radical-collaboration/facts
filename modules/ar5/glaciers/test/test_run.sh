#!/bin/bash

python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --scenario ssp585 --baseyear 2005 --tlm_data 1
python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --gmip 2
python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --chunksize 662 --locationfile $LFILE

