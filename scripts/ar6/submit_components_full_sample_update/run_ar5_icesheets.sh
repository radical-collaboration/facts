#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=ar5ice
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_ar5ice.txt
#SBATCH --error=error_ar5ice.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/ar5icesheets

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "AR5 icesheets - SSP119"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp119 --scenario ssp119 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp119
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp119 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp119 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - SSP126"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp126 --scenario ssp126 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp126
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp126 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - SSP245"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp245 --scenario ssp245 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp245
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp245 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - SSP370"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp370 --scenario ssp370 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp370
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp370 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp370 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - SSP585"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585 --scenario ssp585 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-ssp585
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-ssp585 --chunksize 662 --locationfile ${LFILE}


echo "AR5 icesheets - tlim1.5win0.25"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim1.5win0.25 --scenario tlim1.5win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim1.5win0.25
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim1.5win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim1.5win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - tlim2.0win0.25"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim2.0win0.25 --scenario tlim2.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim2.0win0.25
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim2.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - tlim3.0win0.25"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim3.0win0.25 --scenario tlim3.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim3.0win0.25
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim3.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim3.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - tlim4.0win0.25"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim4.0win0.25 --scenario tlim4.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim4.0win0.25
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim4.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim4.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 icesheets - tlim5.0win0.25"
srun python ar5_preprocess_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim5.0win0.25 --scenario tlim5.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_icesheets.py --pipeline_id icesheets-ar5-icesheets-tlim5.0win0.25
srun python ar5_project_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim5.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ar5_postprocess_icesheets_dev.py --pipeline_id icesheets-ar5-icesheets-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
