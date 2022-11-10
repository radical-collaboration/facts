#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=sej
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_sej.txt
#SBATCH --error=error_sej.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/bambericesheets

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "SEJ icesheets - SSP126"
srun python bamber19_preprocess_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp126 --scenario rcp26 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python bamber19_fit_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp126
srun python bamber19_project_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp126 --nsamps 20000 --seed 4321
srun python bamber19_postprocess_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "SEJ icesheets - SSP245"
srun python bamber19_preprocess_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp245 --scenario rcp45 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python bamber19_fit_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp245
srun python bamber19_project_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp245 --nsamps 20000 --seed 4321
srun python bamber19_postprocess_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "SEJ icesheets - SSP585"
srun python bamber19_preprocess_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp585 --scenario rcp85 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python bamber19_fit_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp585
srun python bamber19_project_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp585 --nsamps 20000 --seed 4321
srun python bamber19_postprocess_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-ssp585 --chunksize 662 --locationfile ${LFILE}

echo "SEJ icesheets - tlim2.0win0.25"
srun python bamber19_preprocess_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim2.0win0.25 --scenario tlim2.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python bamber19_fit_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim2.0win0.25
srun python bamber19_project_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim2.0win0.25 --nsamps 20000 --seed 4321
srun python bamber19_postprocess_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "SEJ icesheets - tlim5.0win0.25"
srun python bamber19_preprocess_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim5.0win0.25 --scenario tlim5.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python bamber19_fit_icesheets.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim5.0win0.25
srun python bamber19_project_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim5.0win0.25 --nsamps 20000 --seed 4321
srun python bamber19_postprocess_icesheets_dev.py --pipeline_id icesheets-ipccar6-bambericesheet-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
