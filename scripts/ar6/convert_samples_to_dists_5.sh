#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=conv_dists5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_convert_samples_to_dists5.txt
#SBATCH --error=error_convert_samples_to_dists5.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"
TEMPDIR="/scratch/ggg46/temp_convert_dists/temp5"

cd ${SCRIPTDIR}


echo "Generating distribution files for the regional components"
srun python convert_samples_to_dist.py --indir ../regional/full_sample_components/ --tempdir ${TEMPDIR} --outdir ../regional/dist_components/ --chunksize 662

#echo "Generating distribution files for the regional components rates"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_components_rates/ --tempdir ${TEMPDIR} --outdir ../regional/dist_components_rates/ --chunksize 662
