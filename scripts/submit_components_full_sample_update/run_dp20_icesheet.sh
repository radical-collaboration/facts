#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=mici
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_mici.txt
#SBATCH --error=error_mici.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/dp20icesheet

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "MICI - SSP126"
srun python dp20_preprocess_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp126 --scenario rcp26 --baseyear 2005
srun python dp20_fit_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp126
srun python dp20_project_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp126 --nsamps 20000 --seed 4321 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python dp20_postprocess_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "MICI - SSP245"
srun python dp20_preprocess_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp245 --scenario rcp45 --baseyear 2005
srun python dp20_fit_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp245
srun python dp20_project_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp245 --nsamps 20000 --seed 4321 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python dp20_postprocess_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "MICI - SSP585"
srun python dp20_preprocess_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp585 --scenario rcp85 --baseyear 2005
srun python dp20_fit_icesheet.py --pipeline_id icesheets-dp20-icesheet-ssp585
srun python dp20_project_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp585 --nsamps 20000 --seed 4321 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python dp20_postprocess_icesheet_dev.py --pipeline_id icesheets-dp20-icesheet-ssp585 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
