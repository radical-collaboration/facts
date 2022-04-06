#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=larmip
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_larmip.txt
#SBATCH --error=error_larmip.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/larmipicesheet

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "LARMIP - SSP119"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp119 --scenario ssp119 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp119
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp119 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp119 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - SSP126"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp126 --scenario ssp126 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp126
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp126 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - SSP245"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp245 --scenario ssp245 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp245
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp245 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - SSP370"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --scenario ssp370 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - SSP585"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp585 --scenario ssp585 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp585
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp585 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp585 --chunksize 662 --locationfile ${LFILE}


echo "LARMIP - tlim1.5win0.25"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim1.5win0.25 --scenario tlim1.5win0.25 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim1.5win0.25
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim1.5win0.25 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim1.5win0.25 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - tlim2.0win0.25"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim2.0win0.25 --scenario tlim2.0win0.25 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim2.0win0.25
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim2.0win0.25 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - tlim3.0win0.25"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim3.0win0.25 --scenario tlim3.0win0.25 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim3.0win0.25
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim3.0win0.25 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim3.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - tlim4.0win0.25"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim4.0win0.25 --scenario tlim4.0win0.25 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim4.0win0.25
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim4.0win0.25 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim4.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "LARMIP - tlim5.0win0.25"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim5.0win0.25 --scenario tlim5.0win0.25 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim5.0win0.25
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim5.0win0.25 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
