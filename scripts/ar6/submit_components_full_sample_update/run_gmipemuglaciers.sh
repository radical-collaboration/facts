#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=gmipemu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_gmip.txt
#SBATCH --error=error_gmip.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/gmipemuglaciers

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


# Original seed - 1234
# Seed to correlate input temperatures - 4321

echo "GMIPemu glaciers - SSP119"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp119 --scenario ssp119 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp119
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp119 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp119 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - SSP126"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp126 --scenario ssp126 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp126
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp126 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - SSP245"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp245 --scenario ssp245 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp245
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp245 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - SSP370"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp370 --scenario ssp370 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp370
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp370 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp370 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - SSP585"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp585 --scenario ssp585 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp585
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp585 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-ssp585 --chunksize 662 --locationfile ${LFILE}


echo "GMIPemu glaciers - tlim1.5win0.25"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25 --scenario tlim1.5win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - tlim2.0win0.25"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25 --scenario tlim2.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - tlim3.0win0.25"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25 --scenario tlim3.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - tlim4.0win0.25"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25 --scenario tlim4.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "GMIPemu glaciers - tlim5.0win0.25"
srun python ipccar6_preprocess_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25 --scenario tlim5.0win0.25 --pyear_start 2020 --pyear_end 2100 --pyear_step 10 --baseyear 2005 --model_driver FAIR
srun python ipccar6_fit_gmipemuglaciers.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25
srun python ipccar6_project_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25 --nsamps 20000 --seed 4321
srun python ipccar6_postprocess_gmipemuglaciers_dev.py --pipeline_id glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
