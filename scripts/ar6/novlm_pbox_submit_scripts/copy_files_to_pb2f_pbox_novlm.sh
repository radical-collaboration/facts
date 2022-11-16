#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/regional_novlm/dist_components"
WFDIR="/projects/kopp/ar6/regional_novlm/dist_workflows"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/regional_novlm/pboxes/pb_2f"

# Go to the pbox directory
cd ${PBDIR}

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_localsl.nc ./ssp126/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp126_GIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_GIS_localsl.nc --outfile ./ssp126/icesheets-pb2f-icesheets-ssp126_GIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp126_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp126_AIS_localsl.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp126_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_AIS_localsl.nc --outfile ./ssp126/icesheets-pb2f-icesheets-ssp126_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp126/total-workflow.nc ${WFDIR}/wf_2f/ssp126/total-workflow.nc ${WFDIR}/wf_3f/ssp126/total-workflow.nc ${WFDIR}/wf_4/ssp126/total-workflow.nc --outfile ./ssp126/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./ssp245/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp245_GIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_GIS_localsl.nc --outfile ./ssp245/icesheets-pb2f-icesheets-ssp245_GIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp245_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp245_AIS_localsl.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp245_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_AIS_localsl.nc --outfile ./ssp245/icesheets-pb2f-icesheets-ssp245_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp245/total-workflow.nc ${WFDIR}/wf_2f/ssp245/total-workflow.nc ${WFDIR}/wf_3f/ssp245/total-workflow.nc ${WFDIR}/wf_4/ssp245/total-workflow.nc --outfile ./ssp245/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_localsl.nc ./ssp585/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp585_GIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_GIS_localsl.nc --outfile ./ssp585/icesheets-pb2f-icesheets-ssp585_GIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp585_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp585_AIS_localsl.nc ${PROJDIR}/icesheets-dp20-icesheet-ssp585_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_AIS_localsl.nc --outfile ./ssp585/icesheets-pb2f-icesheets-ssp585_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp585/total-workflow.nc ${WFDIR}/wf_2f/ssp585/total-workflow.nc ${WFDIR}/wf_3f/ssp585/total-workflow.nc ${WFDIR}/wf_4/ssp585/total-workflow.nc --outfile ./ssp585/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
