#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=ismipemu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_ismipemu.txt
#SBATCH --error=error_ismipemu.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/ismipemuicesheet

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "ISMIPemu icesheets - SSP119"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp119 --scenario ssp119 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp119
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp119 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp119 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - SSP126"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp126 --scenario ssp126 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp126
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp126 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - SSP245"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp245 --scenario ssp245 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp245
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp245 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - SSP370"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp370 --scenario ssp370 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp370
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp370 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp370 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - SSP585"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --scenario ssp585 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-ssp585 --chunksize 662 --locationfile ${LFILE}


echo "ISMIPemu icesheets - tlim1.5win0.25"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25 --scenario tlim1.5win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - tlim2.0win0.25"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25 --scenario tlim2.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - tlim3.0win0.25"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25 --scenario tlim3.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - tlim4.0win0.25"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25 --scenario tlim4.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "ISMIPemu icesheets - tlim5.0win0.25"
srun python ipccar6_preprocess_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25 --scenario tlim5.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_ismipemuicesheet.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25
srun python ipccar6_project_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_ismipemuicesheet_dev.py --pipeline_id icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
