#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=vlm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=96:00:00
#SBATCH --output=out_vlm.txt
#SBATCH --error=error_vlm.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/verticallandmotion

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"


echo "VLM - Running"
srun python kopp14_preprocess_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
srun python kopp14_fit_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
srun python kopp14_project_verticallandmotion.py --pipeline_id verticallandmotion-kopp14-verticallandmotion
srun python kopp14_postprocess_verticallandmotion_dev.py --pipeline_id verticallandmotion-kopp14-verticallandmotion --nsamps 20000 --seed 7331 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --chunksize 662 --locationfile ${LFILE}


mv *localsl* ../output_local

rm *.pkl
