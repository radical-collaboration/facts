#!/bin/bash

python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --scenario ssp5 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585
python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --chunksize 6619 --locationfile $LFILE

