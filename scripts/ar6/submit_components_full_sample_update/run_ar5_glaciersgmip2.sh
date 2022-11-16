#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=ar5glac
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_ar5glac.txt
#SBATCH --error=error_ar5glac.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/ar5glaciers

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


# Original seed - 1234
# Seed to correlate input temperatures - 4321


echo "AR5 glaciers - SSP119"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp119 --scenario ssp119 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp119 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp119 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp119 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - SSP126"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp126 --scenario ssp126 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp126 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp126 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp126 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - SSP245"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp245 --scenario ssp245 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp245 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp245 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp245 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - SSP370"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp370 --scenario ssp370 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp370 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp370 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp370 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - SSP585"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --scenario ssp585 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-ssp585 --chunksize 662 --locationfile ${LFILE}


echo "AR5 glaciers - tlim1.5win0.25"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim1.5win0.25 --scenario tlim1.5win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim1.5win0.25 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim1.5win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim1.5win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - tlim2.0win0.25"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim2.0win0.25 --scenario tlim2.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim2.0win0.25 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim2.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim2.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - tlim3.0win0.25"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim3.0win0.25 --scenario tlim3.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim3.0win0.25 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim3.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim3.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - tlim4.0win0.25"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim4.0win0.25 --scenario tlim4.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim4.0win0.25 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim4.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim4.0win0.25 --chunksize 662 --locationfile ${LFILE}

echo "AR5 glaciers - tlim5.0win0.25"
srun python ar5_preprocess_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim5.0win0.25 --scenario tlim5.0win0.25 --baseyear 2005 --tlm_data 1
srun python ar5_fit_glaciers.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim5.0win0.25 --gmip 2
srun python ar5_project_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim5.0win0.25 --nsamps 20000 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --seed 4321
srun python ar5_postprocess_glaciers_dev.py --pipeline_id glaciers-ar5-glaciersgmip2-tlim5.0win0.25 --chunksize 662 --locationfile ${LFILE}


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
