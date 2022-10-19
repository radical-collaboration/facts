#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=conv_dists0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_convert_samples_to_dists0.txt
#SBATCH --error=error_convert_samples_to_dists0.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"
TEMPDIR="/scratch/ggg46/temp_convert_dists/temp0"

cd ${SCRIPTDIR}


#echo "Generating distribution files for regional workflows"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows/wf_1e --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows/wf_1e --chunksize 662

#echo "Generating distribution files for regional workflows rates"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows_rates/wf_1e --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows_rates/wf_1e --chunksize 662


echo "Generating distribution files for regional no-vlm workflows"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows/wf_1e --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows/wf_1e --chunksize 662

echo "Generating distribution files for regional no-vlm workflows rates"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows_rates/wf_1e --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows_rates/wf_1e --chunksize 662
