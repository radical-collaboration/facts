#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=od
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_od.txt
#SBATCH --error=error_od.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/oceandynamics

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


# Original seed - 1463
# Seed to correlate input temperatures - 4321

echo "OD - SSP119"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp119 --scenario ssp119 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp119
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp119 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp119 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - SSP126"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp126 --scenario ssp126 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp126
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp126 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp126 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - SSP245"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp245 --scenario ssp245 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp245
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp245 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp245 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - SSP370"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp370 --scenario ssp370 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp370
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp370 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp370 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - SSP585"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --scenario ssp585 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-ssp585 --chunksize 662 --nsamps 20000 --seed 4321


echo "OD - tlim1.5win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim1.5win0.25 --scenario tlim1.5win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim1.5win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim1.5win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim1.5win0.25 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - tlim2.0win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim2.0win0.25 --scenario tlim2.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim2.0win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim2.0win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim2.0win0.25 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - tlim3.0win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim3.0win0.25 --scenario tlim3.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim3.0win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim3.0win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim3.0win0.25 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - tlim4.0win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --scenario tlim4.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --chunksize 662 --nsamps 20000 --seed 4321

echo "OD - tlim5.0win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim5.0win0.25 --scenario tlim5.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim5.0win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim5.0win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim5.0win0.25 --chunksize 662 --nsamps 20000 --seed 4321


mv *globalsl* ../output_global
mv *localsl* ../output_local

rm *.pkl
