#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/regional_novlm/dist_components"
WFDIR="/projects/kopp/ar6/regional_novlm/dist_workflows"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/regional_novlm/pboxes/pb_1f"

# Go to the pbox directory
cd ${PBDIR}

# Remove the old files
mkdir -p ssp119
rm ./ssp119/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp119_localsl.nc ./ssp119/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp119_localsl.nc ./ssp119/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp119_localsl.nc ./ssp119/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp119_GIS_localsl.nc ./ssp119/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp119_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp119_AIS_localsl.nc --outfile ./ssp119/icesheets-pb1f-icesheets-ssp119_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp119/total-workflow.nc ${WFDIR}/wf_2f/ssp119/total-workflow.nc --outfile ./ssp119/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp126_GIS_localsl.nc ./ssp126/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp126_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp126_AIS_localsl.nc --outfile ./ssp126/icesheets-pb1f-icesheets-ssp126_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp126/total-workflow.nc ${WFDIR}/wf_2f/ssp126/total-workflow.nc --outfile ./ssp126/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp245_GIS_localsl.nc ./ssp245/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp245_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp245_AIS_localsl.nc --outfile ./ssp245/icesheets-pb1f-icesheets-ssp245_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp245/total-workflow.nc ${WFDIR}/wf_2f/ssp245/total-workflow.nc --outfile ./ssp245/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp370
rm ./ssp370/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp370_localsl.nc ./ssp370/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp370_localsl.nc ./ssp370/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp370_localsl.nc ./ssp370/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp370_GIS_localsl.nc ./ssp370/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp370_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp370_AIS_localsl.nc --outfile ./ssp370/icesheets-pb1f-icesheets-ssp370_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp370/total-workflow.nc ${WFDIR}/wf_2f/ssp370/total-workflow.nc --outfile ./ssp370/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp585_GIS_localsl.nc ./ssp585/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp585_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp585_AIS_localsl.nc --outfile ./ssp585/icesheets-pb1f-icesheets-ssp585_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp585/total-workflow.nc ${WFDIR}/wf_2f/ssp585/total-workflow.nc --outfile ./ssp585/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10



# Remove the old files
mkdir -p tlim1.5win0.25
rm ./tlim1.5win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim1.5win0.25_localsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim1.5win0.25_localsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim1.5win0.25_GIS_localsl.nc ./tlim1.5win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim1.5win0.25_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim1.5win0.25_AIS_localsl.nc --outfile ./tlim1.5win0.25/icesheets-pb1f-icesheets-tlim1.5win0.25_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim1.5win0.25/total-workflow.nc ${WFDIR}/wf_2f/tlim1.5win0.25/total-workflow.nc --outfile ./tlim1.5win0.25/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim2.0win0.25
rm ./tlim2.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim2.0win0.25_localsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim2.0win0.25_localsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim2.0win0.25_GIS_localsl.nc ./tlim2.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim2.0win0.25_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim2.0win0.25_AIS_localsl.nc --outfile ./tlim2.0win0.25/icesheets-pb1f-icesheets-tlim2.0win0.25_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim2.0win0.25/total-workflow.nc ${WFDIR}/wf_2f/tlim2.0win0.25/total-workflow.nc --outfile ./tlim2.0win0.25/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim3.0win0.25
rm ./tlim3.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim3.0win0.25_localsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim3.0win0.25_localsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim3.0win0.25_GIS_localsl.nc ./tlim3.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim3.0win0.25_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim3.0win0.25_AIS_localsl.nc --outfile ./tlim3.0win0.25/icesheets-pb1f-icesheets-tlim3.0win0.25_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim3.0win0.25/total-workflow.nc ${WFDIR}/wf_2f/tlim3.0win0.25/total-workflow.nc --outfile ./tlim3.0win0.25/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim4.0win0.25
rm ./tlim4.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim4.0win0.25_localsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim4.0win0.25_localsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim4.0win0.25_GIS_localsl.nc ./tlim4.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim4.0win0.25_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim4.0win0.25_AIS_localsl.nc --outfile ./tlim4.0win0.25/icesheets-pb1f-icesheets-tlim4.0win0.25_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim4.0win0.25/total-workflow.nc ${WFDIR}/wf_2f/tlim4.0win0.25/total-workflow.nc --outfile ./tlim4.0win0.25/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim5.0win0.25
rm ./tlim5.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim5.0win0.25_localsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim5.0win0.25_localsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim5.0win0.25_GIS_localsl.nc ./tlim5.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim5.0win0.25_AIS_localsl.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim5.0win0.25_AIS_localsl.nc --outfile ./tlim5.0win0.25/icesheets-pb1f-icesheets-tlim5.0win0.25_AIS_localsl.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim5.0win0.25/total-workflow.nc ${WFDIR}/wf_2f/tlim5.0win0.25/total-workflow.nc --outfile ./tlim5.0win0.25/total-workflow.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
