#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=fitismip
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_fitismip.txt
#SBATCH --error=error_fitismip.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/fittedismip

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "FittedISMIP - SSP119"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp119 --scenario ssp119 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp119
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp119 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp119 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - SSP126"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp126 --scenario ssp126 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp126
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp126 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp126 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - SSP245"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp245 --scenario ssp245 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp245
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp245 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp245 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - SSP370"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp370 --scenario ssp370 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp370
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp370 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp370 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - SSP585"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp585 --scenario ssp585 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp585
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp585 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-ssp585 --chunksize 6619 --locationfile ${LFILE}


echo "FittedISMIP - tlim1.5win0.25"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim1.5win0.25 --scenario tlim1.5win0.25 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim1.5win0.25
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim1.5win0.25 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim1.5win0.25 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - tlim2.0win0.25"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim2.0win0.25 --scenario tlim2.0win0.25 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim2.0win0.25
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim2.0win0.25 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim2.0win0.25 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - tlim3.0win0.25"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim3.0win0.25 --scenario tlim3.0win0.25 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim3.0win0.25
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim3.0win0.25 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim3.0win0.25 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - tlim4.0win0.25"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim4.0win0.25 --scenario tlim4.0win0.25 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim4.0win0.25
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim4.0win0.25 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim4.0win0.25 --chunksize 6619 --locationfile ${LFILE}

echo "FittedISMIP - tlim5.0win0.25"
srun python FittedISMIP_preprocess_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim5.0win0.25 --scenario tlim5.0win0.25 --tlm_data 1
srun python FittedISMIP_fit_icesheets.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim5.0win0.25
srun python FittedISMIP_project_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim5.0win0.25 --nsamps 20000 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python FittedISMIP_postprocess_icesheets_dev.py --pipeline_id icesheets-FittedISMIP-icesheets-tlim5.0win0.25 --chunksize 6619 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
