#!/bin/bash
#SBATCH --partition=mem
#SBATCH --job-name=tot_wf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=998000
#SBATCH --time=72:00:00
#SBATCH --output=out_total_workflows2.txt
#SBATCH --error=error_total_workflows2.txt
#SBATCH --export=ALL



######## DO NOT RUN IN PARALLEL WITH GLOBAL TOTALING SCRIPT #############


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

# Define directories
TOTDIR="/scratch/ggg46/facts_total_module"
WFDIR="/projects/kopp/ar6/regional/full_sample_workflows"
TEMPDIR="/scratch/ggg46/temp_total_workflow"

# Go to the total module directory
cd ${TOTDIR}

# Total up the workflows
#echo "1e workflow"

#srun cp -r ${WFDIR}/wf_1e/* ${TEMPDIR}

#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp119 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp370 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2100 --chunksize 1324

#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim1.5win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim2.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim3.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim4.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim5.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324

#rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_1e/
#rm ${TEMPDIR}/*/*.nc
#rmdir ${TEMPDIR}/*


#echo "2e workflow"

#srun cp -r ${WFDIR}/wf_2e/* ${TEMPDIR}

#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp119 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp370 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2100 --chunksize 1324

#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim1.5win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim2.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim3.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim4.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim5.0win0.25 --pyear_start 2020 --pyear_end 2100 --chunksize 1324

#rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_2e/
#rm ${TEMPDIR}/*/*.nc
#rmdir ${TEMPDIR}/*


#echo "3e workflow"

#srun cp -r ${WFDIR}/wf_3e/* ${TEMPDIR}

#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2100 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2100 --chunksize 1324

#rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_3e/
#rm ${TEMPDIR}/*/*.nc
#rmdir ${TEMPDIR}/*


echo "1f workflow"

#srun cp -r ${WFDIR}/wf_1f/* ${TEMPDIR}

#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp119 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp370 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim1.5win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim2.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim3.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
#srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim4.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim5.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_1f/
rm ${TEMPDIR}/*/*.nc
rmdir ${TEMPDIR}/*


echo "2f workflow"

srun cp -r ${WFDIR}/wf_2f/* ${TEMPDIR}

srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp119 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp370 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim1.5win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim2.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim3.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim4.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim5.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_2f/
rm ${TEMPDIR}/*/*.nc
rmdir ${TEMPDIR}/*


echo "3f workflow"

srun cp -r ${WFDIR}/wf_3f/* ${TEMPDIR}

srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_3f/
rm ${TEMPDIR}/*/*.nc
rmdir ${TEMPDIR}/*


echo "4 workflow"

srun cp -r ${WFDIR}/wf_4/* ${TEMPDIR}

srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp126 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp245 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/ssp585 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim2.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324
srun python total_workflow_dev.py --directory ${TEMPDIR}/tlim5.0win0.25 --pyear_start 2020 --pyear_end 2300 --chunksize 1324

rsync -avm --include='total-workflow.nc' -f 'hide,! */' ${TEMPDIR}/ ${WFDIR}/wf_4/
rm ${TEMPDIR}/*/*.nc
rmdir ${TEMPDIR}/*


# Remove the leftover total file from the module directory
echo "Cleaning up..."
rm total-workflow.nc
