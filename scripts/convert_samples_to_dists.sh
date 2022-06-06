#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=conv_dists
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_convert_samples_to_dists.txt
#SBATCH --error=error_convert_samples_to_dists.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"

cd ${SCRIPTDIR}

echo "Generating distribution files for global workflows"
srun python convert_samples_to_dist.py --indir ../global/full_sample_workflows/ --outdir ../global/dist_workflows/ --chunksize 6619

echo "Generating distribution files for the global components"
srun python convert_samples_to_dist.py --indir ../global/full_sample_components/ --outdir ../global/dist_components/ --chunksize 6619

echo "Generating distribution files for regional workflows"
srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows/ --outdir ../regional/dist_workflows/ --chunksize 6619

echo "Generating distribution files for the regional components"
srun python convert_samples_to_dist.py --indir ../regional/full_sample_components/ --outdir ../regional/dist_components/ --chunksize 6619

#echo "Generating distribution files for global workflows rates"
#srun python convert_samples_to_dist.py --indir ../global/full_sample_workflows_rates/ --outdir ../global/dist_workflows_rates/ --chunksize 6619

#echo "Generating distribution files for the global components rates"
#srun python convert_samples_to_dist.py --indir ../global/full_sample_components_rates/ --outdir ../global/dist_components_rates/ --chunksize 6619

#echo "Generating distribution files for regional workflows rates"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows_rates/ --outdir ../regional/dist_workflows_rates/ --chunksize 6619

#echo "Generating distribution files for the regional components rates"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_components_rates/ --outdir ../regional/dist_components_rates/ --chunksize 6619
