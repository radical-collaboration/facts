#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=gen_rates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_generate_rate_files.txt
#SBATCH --error=error_generate_rate_files.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"

cd ${SCRIPTDIR}

echo "Generating rate files for global workflows"
srun python generate_rate_files.py --indir ../global/full_sample_workflows/ --outdir ../global/full_sample_workflows_rates/ --chunksize 662

echo "Generating rate files for the global components"
srun python generate_rate_files.py --indir ../global/full_sample_components/ --outdir ../global/full_sample_components_rates/ --chunksize 662

echo "Generating rate files for regional workflows"
srun python generate_rate_files.py --indir ../regional/full_sample_workflows/ --outdir ../regional/full_sample_workflows_rates/ --chunksize 662

echo "Generating rate files for the regional components"
srun python generate_rate_files.py --indir ../regional/full_sample_components/ --outdir ../regional/full_sample_components_rates/ --chunksize 662
