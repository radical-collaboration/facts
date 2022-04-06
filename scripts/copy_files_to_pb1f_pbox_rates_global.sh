#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/dist_components_rates"
WFDIR="/projects/kopp/ar6/global/dist_workflows_rates"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/global/pboxes_rates/pb_1f"

# Go to the pbox directory
cd ${PBDIR}

# Remove the old files
mkdir -p ssp119
rm ./ssp119/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp119_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp119_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp119_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp119_GIS_globalsl_rates.nc ./ssp119/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp119_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp119_TOT_globalsl_rates.nc --outfile ./ssp119/icesheets-pb1f-icesheets-ssp119_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp119/total-workflow_rates.nc ${WFDIR}/wf_2f/ssp119/total-workflow_rates.nc --outfile ./ssp119/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp126_GIS_globalsl_rates.nc ./ssp126/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp126_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp126_TOT_globalsl_rates.nc --outfile ./ssp126/icesheets-pb1f-icesheets-ssp126_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp126/total-workflow_rates.nc ${WFDIR}/wf_2f/ssp126/total-workflow_rates.nc --outfile ./ssp126/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp245_GIS_globalsl_rates.nc ./ssp245/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp245_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp245_TOT_globalsl_rates.nc --outfile ./ssp245/icesheets-pb1f-icesheets-ssp245_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp245/total-workflow_rates.nc ${WFDIR}/wf_2f/ssp245/total-workflow_rates.nc --outfile ./ssp245/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp370
rm ./ssp370/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp370_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp370_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp370_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp370_GIS_globalsl_rates.nc ./ssp370/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp370_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp370_TOT_globalsl_rates.nc --outfile ./ssp370/icesheets-pb1f-icesheets-ssp370_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp370/total-workflow_rates.nc ${WFDIR}/wf_2f/ssp370/total-workflow_rates.nc --outfile ./ssp370/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp585_GIS_globalsl_rates.nc ./ssp585/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-ssp585_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp585_TOT_globalsl_rates.nc --outfile ./ssp585/icesheets-pb1f-icesheets-ssp585_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/ssp585/total-workflow_rates.nc ${WFDIR}/wf_2f/ssp585/total-workflow_rates.nc --outfile ./ssp585/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10



# Remove the old files
mkdir -p tlim1.5win0.25
rm ./tlim1.5win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim1.5win0.25_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim1.5win0.25_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim1.5win0.25_GIS_globalsl_rates.nc ./tlim1.5win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim1.5win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim1.5win0.25_TOT_globalsl_rates.nc --outfile ./tlim1.5win0.25/icesheets-pb1f-icesheets-tlim1.5win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim1.5win0.25/total-workflow_rates.nc ${WFDIR}/wf_2f/tlim1.5win0.25/total-workflow_rates.nc --outfile ./tlim1.5win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim2.0win0.25
rm ./tlim2.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim2.0win0.25_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim2.0win0.25_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim2.0win0.25_GIS_globalsl_rates.nc ./tlim2.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim2.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim2.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim2.0win0.25/icesheets-pb1f-icesheets-tlim2.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim2.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2f/tlim2.0win0.25/total-workflow_rates.nc --outfile ./tlim2.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim3.0win0.25
rm ./tlim3.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim3.0win0.25_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim3.0win0.25_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim3.0win0.25_GIS_globalsl_rates.nc ./tlim3.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim3.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim3.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim3.0win0.25/icesheets-pb1f-icesheets-tlim3.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim3.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2f/tlim3.0win0.25/total-workflow_rates.nc --outfile ./tlim3.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim4.0win0.25
rm ./tlim4.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim4.0win0.25_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim4.0win0.25_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim4.0win0.25_GIS_globalsl_rates.nc ./tlim4.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim4.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim4.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim4.0win0.25/icesheets-pb1f-icesheets-tlim4.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim4.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2f/tlim4.0win0.25/total-workflow_rates.nc --outfile ./tlim4.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10

# Remove the old files
mkdir -p tlim5.0win0.25
rm ./tlim5.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim5.0win0.25_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim5.0win0.25_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-tlim5.0win0.25_GIS_globalsl_rates.nc ./tlim5.0win0.25/

# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ar5-icesheets-tlim5.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim5.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim5.0win0.25/icesheets-pb1f-icesheets-tlim5.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1f/tlim5.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2f/tlim5.0win0.25/total-workflow_rates.nc --outfile ./tlim5.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2300 --pyear_step 10
