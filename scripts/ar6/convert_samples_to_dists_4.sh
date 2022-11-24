#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=conv_dists4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=192000
#SBATCH --time=144:00:00
#SBATCH --output=out_convert_samples_to_dists4.txt
#SBATCH --error=error_convert_samples_to_dists4.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directories
SCRIPTDIR="/projects/kopp/ar6/scripts"
TEMPDIR="/scratch/ggg46/temp_convert_dists/temp4"

cd ${SCRIPTDIR}


#echo "Generating distribution files for regional workflow 3e"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows/wf_3e --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows/wf_3e --chunksize 662

#echo "Generating distribution files for regional workflow 3f"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows/wf_3f --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows/wf_3f --chunksize 662

#echo "Generating distribution files for regional workflow 4"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows/wf_4 --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows/wf_4 --chunksize 662

#echo "Generating distribution files for regional workflows rates 3e"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows_rates/wf_3e --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows_rates/wf_3e --chunksize 662

#echo "Generating distribution files for regional workflows rates 3f"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows_rates/wf_3f --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows_rates/wf_3f --chunksize 662

#echo "Generating distribution files for regional workflows rates 4"
#srun python convert_samples_to_dist.py --indir ../regional/full_sample_workflows_rates/wf_4 --tempdir ${TEMPDIR} --outdir ../regional/dist_workflows_rates/wf_4 --chunksize 662


echo "Generating distribution files for regional no-vlm workflow 3e"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows/wf_3e --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows/wf_3e --chunksize 662

echo "Generating distribution files for regional no-vlm workflow 3f"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows/wf_3f --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows/wf_3f --chunksize 662

echo "Generating distribution files for regional no-vlm workflow 4"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows/wf_4 --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows/wf_4 --chunksize 662

echo "Generating distribution files for regional no-vlm workflows rates 3e"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows_rates/wf_3e --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows_rates/wf_3e --chunksize 662

echo "Generating distribution files for regional no-vlm workflows rates 3f"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows_rates/wf_3f --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows_rates/wf_3f --chunksize 662

echo "Generating distribution files for regional no-vlm workflows rates 4"
srun python convert_samples_to_dist.py --indir ../regional_novlm/full_sample_workflows_rates/wf_4 --tempdir ${TEMPDIR} --outdir ../regional_novlm/dist_workflows_rates/wf_4 --chunksize 662
