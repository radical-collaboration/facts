#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=gen_rates2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_generate_rate_files2.txt
#SBATCH --error=error_generate_rate_files2.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"
TEMPDIR="/scratch/ggg46/temp_convert_rates/temp2"

cd ${SCRIPTDIR}

echo "Generating rate files for regional workflow"
srun python generate_rate_files.py --indir ../regional/full_sample_workflows/wf_2e --tempdir ${TEMPDIR} --outdir ../regional/full_sample_workflows_rates/wf_2e --chunksize 6619
