#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=larmip3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=8:00:00
#SBATCH --output=out_larmip3.txt
#SBATCH --error=error_larmip3.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/larmipicesheet

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "LARMIP - SSP370"
srun python ipccar6_preprocess_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --scenario ssp370 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --baseyear 2005 --tlm_data 1
srun python ipccar6_fit_larmipicesheet.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370
srun python ipccar6_project_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --nsamps 20000 --seed 4321 --crateyear_start 2080 --crateyear_end 2100
srun python ipccar6_postprocess_larmipicesheet_dev.py --pipeline_id icesheets-ipccar6-larmipicesheet-ssp370 --chunksize 662 --locationfile ${LFILE}


mv *ssp370*globalsl* ../output_global
mv *ssp370*localsl* ../output_local

rm *ssp370*.pkl
