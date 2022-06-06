#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=lws
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_lws.txt
#SBATCH --error=error_lws.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/landwaterstorage

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "LWS - SSP119"
srun python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp119 --scenario ssp1 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
srun python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp119
srun python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp119 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
srun python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp119 --chunksize 6619 --locationfile ${LFILE}

echo "LWS - SSP126"
srun python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp126 --scenario ssp1 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
srun python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp126
srun python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp126 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
srun python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp126 --chunksize 6619 --locationfile ${LFILE}

echo "LWS - SSP245"
srun python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp245 --scenario ssp2 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
srun python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp245
srun python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp245 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
srun python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp245 --chunksize 6619 --locationfile ${LFILE}

echo "LWS - SSP370"
srun python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp370 --scenario ssp3 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
srun python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp370
srun python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp370 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
srun python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp370 --chunksize 6619 --locationfile ${LFILE}

echo "LWS - SSP585"
srun python ssp_preprocess_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --scenario ssp5 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005
srun python ssp_fit_landwaterstorage.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585
srun python ssp_project_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --nsamps 20000 --seed 1337 --dcrate_lo -0.4
srun python ssp_postprocess_landwaterstorage_dev.py --pipeline_id landwaterstorage-ssp-landwaterstorage-ssp585 --chunksize 6619 --locationfile ${LFILE}

mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
