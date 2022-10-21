#!/bin/bash

python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --scenario ssp585 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile $LFILE
python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585
python tlm_project_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --nsamps 20000 --seed 4321 --tlm 1
python tlm_postprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --chunksize 662 --nsamps 20000 --seed 4321

