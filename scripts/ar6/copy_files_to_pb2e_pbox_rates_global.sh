#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/dist_components_rates"
WFDIR="/projects/kopp/ar6/global/dist_workflows_rates"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/global/pboxes_rates/pb_2e"

# Go to the pbox directory
cd ${PBDIR}

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl_rates.nc ./ssp126/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp126_globalsl_rates.nc ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_globalsl_rates.nc --outfile ./ssp126/glaciers-pb2e-glaciers-ssp126_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_GIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_GIS_globalsl_rates.nc --outfile ./ssp126/icesheets-pb2e-icesheets-ssp126_GIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp126_TOT_globalsl_rates.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp126_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_AIS_globalsl_rates.nc --outfile ./ssp126/icesheets-pb2e-icesheets-ssp126_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp126/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp126/total-workflow_rates.nc ${WFDIR}/wf_3e/ssp126/total-workflow_rates.nc ${WFDIR}/wf_4/ssp126/total-workflow_rates.nc --outfile ./ssp126/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./ssp245/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp245_globalsl_rates.nc ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_globalsl_rates.nc --outfile ./ssp245/glaciers-pb2e-glaciers-ssp245_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_GIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_GIS_globalsl_rates.nc --outfile ./ssp245/icesheets-pb2e-icesheets-ssp245_GIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp245_TOT_globalsl_rates.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp245_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_AIS_globalsl_rates.nc --outfile ./ssp245/icesheets-pb2e-icesheets-ssp245_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp245/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp245/total-workflow_rates.nc ${WFDIR}/wf_3e/ssp245/total-workflow_rates.nc ${WFDIR}/wf_4/ssp245/total-workflow_rates.nc --outfile ./ssp245/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl_rates.nc ./ssp585/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp585_globalsl_rates.nc ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_globalsl_rates.nc --outfile ./ssp585/glaciers-pb2e-glaciers-ssp585_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_GIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_GIS_globalsl_rates.nc --outfile ./ssp585/icesheets-pb2e-icesheets-ssp585_GIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp585_TOT_globalsl_rates.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp585_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_AIS_globalsl_rates.nc --outfile ./ssp585/icesheets-pb2e-icesheets-ssp585_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp585/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp585/total-workflow_rates.nc ${WFDIR}/wf_3e/ssp585/total-workflow_rates.nc ${WFDIR}/wf_4/ssp585/total-workflow_rates.nc --outfile ./ssp585/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
