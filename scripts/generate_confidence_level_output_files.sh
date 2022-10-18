#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=confid_levs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=out_generate_confidence_level_output_files.txt
#SBATCH --error=error_generate_confidence_level_output_files.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories --------------------------------
MAINDIR="/projects/kopp/ar6"
SCRIPTDIR="${MAINDIR}/scripts"

# Go to the script directory
cd ${SCRIPTDIR}

# Generate the global confidence level files
#echo "Processing global confidence level output files for values"
#srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/global/pboxes/ --outdir ${MAINDIR}/global/confidence_output_files/ --chunksize 66200

# Generate the global confidence level files for rates
#echo "Processing global confidence level output files for rates"
#srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/global/pboxes_rates/ --outdir ${MAINDIR}/global/confidence_output_files/ --chunksize 66200

# Generate the regional confidence level files
#echo "Processing regional confidence level output files for values"
#srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/regional/pboxes/ --outdir ${MAINDIR}/regional/confidence_output_files/ --chunksize 66200

# Generate the regional confidence level files for rates
#echo "Processing regional confidence level output files for rates"
#srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/regional/pboxes_rates/ --outdir ${MAINDIR}/regional/confidence_output_files/ --chunksize 66200

# Generate the regional novlm confidence level files
echo "Processing regional novlm confidence level output files for values"
srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/regional_novlm/pboxes/ --outdir ${MAINDIR}/regional_novlm/confidence_output_files/ --chunksize 66200

# Generate the regional novlm confidence level files for rates
echo "Processing regional novlm confidence level output files for rates"
srun python generate_confidence_output_files.py --pboxdir ${MAINDIR}/regional_novlm/pboxes_rates/ --outdir ${MAINDIR}/regional_novlm/confidence_output_files/ --chunksize 66200
