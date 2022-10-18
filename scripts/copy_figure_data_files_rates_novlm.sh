#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=copy_figs_rates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=144:00:00
#SBATCH --output=out_copy_figure_data_files_rates.txt
#SBATCH --error=error_copy_figure_data_files_rates.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Workflow and projection component directories --------------------------------
FIGDIR="/projects/kopp/ar6/regional_novlm/figure_data_rates"
WFDIR="/projects/kopp/ar6/regional_novlm/dist_workflows_rates"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/regional_novlm/pboxes_rates"

# Go to the pbox directory
cd ${SCRIPTDIR}

# Remove the old figure data files
rm ${FIGDIR}/*/*/*.nc

# Generate figure data by scenario for each workflow and pbox
echo "Processing: SSP119"
mkdir -p ${FIGDIR}/pb_1f/ssp119/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/ssp119/ --outdir ${FIGDIR}/pb_1f/ssp119/
mkdir -p ${FIGDIR}/wf_1f/ssp119/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/ssp119/ --outdir ${FIGDIR}/wf_1f/ssp119/
mkdir -p ${FIGDIR}/wf_2f/ssp119/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/ssp119/ --outdir ${FIGDIR}/wf_2f/ssp119/

mkdir -p ${FIGDIR}/pb_1e/ssp119/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/ssp119/ --outdir ${FIGDIR}/pb_1e/ssp119/
mkdir -p ${FIGDIR}/wf_1e/ssp119/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/ssp119/ --outdir ${FIGDIR}/wf_1e/ssp119/
mkdir -p ${FIGDIR}/wf_2e/ssp119/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/ssp119/ --outdir ${FIGDIR}/wf_2e/ssp119/

echo "Processing: SSP126"
mkdir -p ${FIGDIR}/pb_1f/ssp126/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/ssp126/ --outdir ${FIGDIR}/pb_1f/ssp126/
mkdir -p ${FIGDIR}/pb_2f/ssp126/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2f/ssp126/ --outdir ${FIGDIR}/pb_2f/ssp126/
mkdir -p ${FIGDIR}/wf_1f/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/ssp126/ --outdir ${FIGDIR}/wf_1f/ssp126/
mkdir -p ${FIGDIR}/wf_2f/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/ssp126/ --outdir ${FIGDIR}/wf_2f/ssp126/
mkdir -p ${FIGDIR}/wf_3f/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3f/ssp126/ --outdir ${FIGDIR}/wf_3f/ssp126/

mkdir -p ${FIGDIR}/pb_1e/ssp126/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/ssp126/ --outdir ${FIGDIR}/pb_1e/ssp126/
mkdir -p ${FIGDIR}/pb_2e/ssp126/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2e/ssp126/ --outdir ${FIGDIR}/pb_2e/ssp126/
mkdir -p ${FIGDIR}/wf_1e/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/ssp126/ --outdir ${FIGDIR}/wf_1e/ssp126/
mkdir -p ${FIGDIR}/wf_2e/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/ssp126/ --outdir ${FIGDIR}/wf_2e/ssp126/
mkdir -p ${FIGDIR}/wf_3e/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3e/ssp126/ --outdir ${FIGDIR}/wf_3e/ssp126/

mkdir -p ${FIGDIR}/wf_4/ssp126/
srun python extract_figure_data.py --indir ${WFDIR}/wf_4/ssp126/ --outdir ${FIGDIR}/wf_4/ssp126/

echo "Processing: SSP245"
mkdir -p ${FIGDIR}/pb_1f/ssp245/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/ssp245/ --outdir ${FIGDIR}/pb_1f/ssp245/
mkdir -p ${FIGDIR}/pb_2f/ssp245/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2f/ssp245/ --outdir ${FIGDIR}/pb_2f/ssp245/
mkdir -p ${FIGDIR}/wf_1f/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/ssp245/ --outdir ${FIGDIR}/wf_1f/ssp245/
mkdir -p ${FIGDIR}/wf_2f/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/ssp245/ --outdir ${FIGDIR}/wf_2f/ssp245/
mkdir -p ${FIGDIR}/wf_3f/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3f/ssp245/ --outdir ${FIGDIR}/wf_3f/ssp245/

mkdir -p ${FIGDIR}/pb_1e/ssp245/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/ssp245/ --outdir ${FIGDIR}/pb_1e/ssp245/
mkdir -p ${FIGDIR}/pb_2e/ssp245/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2e/ssp245/ --outdir ${FIGDIR}/pb_2e/ssp245/
mkdir -p ${FIGDIR}/wf_1e/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/ssp245/ --outdir ${FIGDIR}/wf_1e/ssp245/
mkdir -p ${FIGDIR}/wf_2e/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/ssp245/ --outdir ${FIGDIR}/wf_2e/ssp245/
mkdir -p ${FIGDIR}/wf_3e/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3e/ssp245/ --outdir ${FIGDIR}/wf_3e/ssp245/

mkdir -p ${FIGDIR}/wf_4/ssp245/
srun python extract_figure_data.py --indir ${WFDIR}/wf_4/ssp245/ --outdir ${FIGDIR}/wf_4/ssp245/

echo "Processing: SSP370"
mkdir -p ${FIGDIR}/pb_1f/ssp370/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/ssp370/ --outdir ${FIGDIR}/pb_1f/ssp370/
mkdir -p ${FIGDIR}/wf_1f/ssp370/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/ssp370/ --outdir ${FIGDIR}/wf_1f/ssp370/
mkdir -p ${FIGDIR}/wf_2f/ssp370/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/ssp370/ --outdir ${FIGDIR}/wf_2f/ssp370/

mkdir -p ${FIGDIR}/pb_1e/ssp370/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/ssp370/ --outdir ${FIGDIR}/pb_1e/ssp370/
mkdir -p ${FIGDIR}/wf_1e/ssp370/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/ssp370/ --outdir ${FIGDIR}/wf_1e/ssp370/
mkdir -p ${FIGDIR}/wf_2e/ssp370/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/ssp370/ --outdir ${FIGDIR}/wf_2e/ssp370/


echo "Processing: SSP585"
mkdir -p ${FIGDIR}/pb_1f/ssp585/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/ssp585/ --outdir ${FIGDIR}/pb_1f/ssp585/
mkdir -p ${FIGDIR}/pb_2f/ssp585/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2f/ssp585/ --outdir ${FIGDIR}/pb_2f/ssp585/
mkdir -p ${FIGDIR}/wf_1f/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/ssp585/ --outdir ${FIGDIR}/wf_1f/ssp585/
mkdir -p ${FIGDIR}/wf_2f/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/ssp585/ --outdir ${FIGDIR}/wf_2f/ssp585/
mkdir -p ${FIGDIR}/wf_3f/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3f/ssp585/ --outdir ${FIGDIR}/wf_3f/ssp585/

mkdir -p ${FIGDIR}/pb_1e/ssp585/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/ssp585/ --outdir ${FIGDIR}/pb_1e/ssp585/
mkdir -p ${FIGDIR}/pb_2e/ssp585/
srun python extract_figure_data.py --indir ${PBDIR}/pb_2e/ssp585/ --outdir ${FIGDIR}/pb_2e/ssp585/
mkdir -p ${FIGDIR}/wf_1e/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/ssp585/ --outdir ${FIGDIR}/wf_1e/ssp585/
mkdir -p ${FIGDIR}/wf_2e/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/ssp585/ --outdir ${FIGDIR}/wf_2e/ssp585/
mkdir -p ${FIGDIR}/wf_3e/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_3e/ssp585/ --outdir ${FIGDIR}/wf_3e/ssp585/

mkdir -p ${FIGDIR}/wf_4/ssp585/
srun python extract_figure_data.py --indir ${WFDIR}/wf_4/ssp585/ --outdir ${FIGDIR}/wf_4/ssp585/





echo "Processing: tlim1.5win0.25"
mkdir -p ${FIGDIR}/pb_1f/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/tlim1.5win0.25/ --outdir ${FIGDIR}/pb_1f/tlim1.5win0.25/
mkdir -p ${FIGDIR}/wf_1f/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/tlim1.5win0.25/ --outdir ${FIGDIR}/wf_1f/tlim1.5win0.25/
mkdir -p ${FIGDIR}/wf_2f/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/tlim1.5win0.25/ --outdir ${FIGDIR}/wf_2f/tlim1.5win0.25/

mkdir -p ${FIGDIR}/pb_1e/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/tlim1.5win0.25/ --outdir ${FIGDIR}/pb_1e/tlim1.5win0.25/
mkdir -p ${FIGDIR}/wf_1e/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/tlim1.5win0.25/ --outdir ${FIGDIR}/wf_1e/tlim1.5win0.25/
mkdir -p ${FIGDIR}/wf_2e/tlim1.5win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/tlim1.5win0.25/ --outdir ${FIGDIR}/wf_2e/tlim1.5win0.25/

echo "Processing: tlim2.0win0.25"
mkdir -p ${FIGDIR}/pb_1f/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/tlim2.0win0.25/ --outdir ${FIGDIR}/pb_1f/tlim2.0win0.25/
mkdir -p ${FIGDIR}/wf_1f/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/tlim2.0win0.25/ --outdir ${FIGDIR}/wf_1f/tlim2.0win0.25/
mkdir -p ${FIGDIR}/wf_2f/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/tlim2.0win0.25/ --outdir ${FIGDIR}/wf_2f/tlim2.0win0.25/

mkdir -p ${FIGDIR}/pb_1e/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/tlim2.0win0.25/ --outdir ${FIGDIR}/pb_1e/tlim2.0win0.25/
mkdir -p ${FIGDIR}/wf_1e/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/tlim2.0win0.25/ --outdir ${FIGDIR}/wf_1e/tlim2.0win0.25/
mkdir -p ${FIGDIR}/wf_2e/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/tlim2.0win0.25/ --outdir ${FIGDIR}/wf_2e/tlim2.0win0.25/

mkdir -p ${FIGDIR}/wf_4/tlim2.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_4/tlim2.0win0.25/ --outdir ${FIGDIR}/wf_4/tlim2.0win0.25/

echo "Processing: tlim3.0win0.25"
mkdir -p ${FIGDIR}/pb_1f/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/tlim3.0win0.25/ --outdir ${FIGDIR}/pb_1f/tlim3.0win0.25/
mkdir -p ${FIGDIR}/wf_1f/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/tlim3.0win0.25/ --outdir ${FIGDIR}/wf_1f/tlim3.0win0.25/
mkdir -p ${FIGDIR}/wf_2f/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/tlim3.0win0.25/ --outdir ${FIGDIR}/wf_2f/tlim3.0win0.25/

mkdir -p ${FIGDIR}/pb_1e/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/tlim3.0win0.25/ --outdir ${FIGDIR}/pb_1e/tlim3.0win0.25/
mkdir -p ${FIGDIR}/wf_1e/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/tlim3.0win0.25/ --outdir ${FIGDIR}/wf_1e/tlim3.0win0.25/
mkdir -p ${FIGDIR}/wf_2e/tlim3.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/tlim3.0win0.25/ --outdir ${FIGDIR}/wf_2e/tlim3.0win0.25/

echo "Processing: tlim4.0win0.25"
mkdir -p ${FIGDIR}/pb_1f/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/tlim4.0win0.25/ --outdir ${FIGDIR}/pb_1f/tlim4.0win0.25/
mkdir -p ${FIGDIR}/wf_1f/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/tlim4.0win0.25/ --outdir ${FIGDIR}/wf_1f/tlim4.0win0.25/
mkdir -p ${FIGDIR}/wf_2f/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/tlim4.0win0.25/ --outdir ${FIGDIR}/wf_2f/tlim4.0win0.25/

mkdir -p ${FIGDIR}/pb_1e/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/tlim4.0win0.25/ --outdir ${FIGDIR}/pb_1e/tlim4.0win0.25/
mkdir -p ${FIGDIR}/wf_1e/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/tlim4.0win0.25/ --outdir ${FIGDIR}/wf_1e/tlim4.0win0.25/
mkdir -p ${FIGDIR}/wf_2e/tlim4.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/tlim4.0win0.25/ --outdir ${FIGDIR}/wf_2e/tlim4.0win0.25/

echo "Processing: tlim5.0win0.25"
mkdir -p ${FIGDIR}/pb_1f/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1f/tlim5.0win0.25/ --outdir ${FIGDIR}/pb_1f/tlim5.0win0.25/
mkdir -p ${FIGDIR}/wf_1f/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1f/tlim5.0win0.25/ --outdir ${FIGDIR}/wf_1f/tlim5.0win0.25/
mkdir -p ${FIGDIR}/wf_2f/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2f/tlim5.0win0.25/ --outdir ${FIGDIR}/wf_2f/tlim5.0win0.25/

mkdir -p ${FIGDIR}/pb_1e/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${PBDIR}/pb_1e/tlim5.0win0.25/ --outdir ${FIGDIR}/pb_1e/tlim5.0win0.25/
mkdir -p ${FIGDIR}/wf_1e/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_1e/tlim5.0win0.25/ --outdir ${FIGDIR}/wf_1e/tlim5.0win0.25/
mkdir -p ${FIGDIR}/wf_2e/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_2e/tlim5.0win0.25/ --outdir ${FIGDIR}/wf_2e/tlim5.0win0.25/

mkdir -p ${FIGDIR}/wf_4/tlim5.0win0.25/
srun python extract_figure_data.py --indir ${WFDIR}/wf_4/tlim5.0win0.25/ --outdir ${FIGDIR}/wf_4/tlim5.0win0.25/
