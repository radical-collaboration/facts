#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=4remvlm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=192000
#SBATCH --time=144:00:00
#SBATCH --output=out_remove_vlm_from_total4.txt
#SBATCH --error=error_remove_vlm_from_total4.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Directory definitions
MAINDIR="/projects/kopp/ar6"
SCRIPTDIR="${MAINDIR}/scripts"

CHUNKSIZE=6619


# Change to the script directory
cd ${SCRIPTDIR}


# Values =======================================================================
BINDIR="${MAINDIR}/regional/full_sample_workflows/wf_3e"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows/wf_3e"

echo "==== Removing VLM from values ===="
INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}


BINDIR="${MAINDIR}/regional/full_sample_workflows/wf_3f"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows/wf_3f"

INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}


BINDIR="${MAINDIR}/regional/full_sample_workflows/wf_4"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows/wf_4"

INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/tlim2.0win0.25"
OUTDIR="${BOUTDIR}/tlim2.0win0.25"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/tlim5.0win0.25"
OUTDIR="${BOUTDIR}/tlim5.0win0.25"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}


# Rates ========================================================================
BINDIR="${MAINDIR}/regional/full_sample_workflows_rates/wf_3e"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows_rates/wf_3e"

echo "==== Removing VLM from values ===="
INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}


BINDIR="${MAINDIR}/regional/full_sample_workflows_rates/wf_3f"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows_rates/wf_3f"

INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}


BINDIR="${MAINDIR}/regional/full_sample_workflows_rates/wf_4"
BOUTDIR="${MAINDIR}/regional_novlm/full_sample_workflows_rates/wf_4"

INDIR="${BINDIR}/ssp126"
OUTDIR="${BOUTDIR}/ssp126"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp245"
OUTDIR="${BOUTDIR}/ssp245"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/ssp585"
OUTDIR="${BOUTDIR}/ssp585"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/tlim2.0win0.25"
OUTDIR="${BOUTDIR}/tlim2.0win0.25"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}

INDIR="${BINDIR}/tlim5.0win0.25"
OUTDIR="${BOUTDIR}/tlim5.0win0.25"
echo "INDIR: ${INDIR}"
echo "OUTDIR: ${OUTDIR}"
srun python remove_vlm_from_total.py --indir ${INDIR} --outdir ${OUTDIR} --chunksize ${CHUNKSIZE}
