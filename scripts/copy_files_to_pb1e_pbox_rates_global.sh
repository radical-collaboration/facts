#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/dist_components_rates"
WFDIR="/projects/kopp/ar6/global/dist_workflows_rates"
SCRIPTDIR="/projects/kopp/ar6/scripts"
PBDIR="/projects/kopp/ar6/global/pboxes_rates/pb_1e"

# Go to the pbox directory
cd ${PBDIR}

# Remove the old files
mkdir -p ssp119
rm ./ssp119/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp119_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp119_GIS_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp119_globalsl_rates.nc ./ssp119/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp119_globalsl_rates.nc ./ssp119/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp119_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp119_TOT_globalsl_rates.nc --outfile ./ssp119/icesheets-pb1e-icesheets-ssp119_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp119/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp119/total-workflow_rates.nc --outfile ./ssp119/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_GIS_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl_rates.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl_rates.nc ./ssp126/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp126_TOT_globalsl_rates.nc --outfile ./ssp126/icesheets-pb1e-icesheets-ssp126_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp126/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp126/total-workflow_rates.nc --outfile ./ssp126/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_GIS_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl_rates.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./ssp245/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp245_TOT_globalsl_rates.nc --outfile ./ssp245/icesheets-pb1e-icesheets-ssp245_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp245/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp245/total-workflow_rates.nc --outfile ./ssp245/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p ssp370
rm ./ssp370/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp370_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp370_GIS_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp370_globalsl_rates.nc ./ssp370/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp370_globalsl_rates.nc ./ssp370/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp370_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp370_TOT_globalsl_rates.nc --outfile ./ssp370/icesheets-pb1e-icesheets-ssp370_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp370/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp370/total-workflow_rates.nc --outfile ./ssp370/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_GIS_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl_rates.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl_rates.nc ./ssp585/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-ssp585_TOT_globalsl_rates.nc --outfile ./ssp585/icesheets-pb1e-icesheets-ssp585_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/ssp585/total-workflow_rates.nc ${WFDIR}/wf_2e/ssp585/total-workflow_rates.nc --outfile ./ssp585/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10




# Remove the old files
mkdir -p tlim1.5win0.25
rm ./tlim1.5win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25_GIS_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim1.5win0.25_globalsl_rates.nc ./tlim1.5win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim1.5win0.25/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim1.5win0.25_TOT_globalsl_rates.nc --outfile ./tlim1.5win0.25/icesheets-pb1e-icesheets-tlim1.5win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/tlim1.5win0.25/total-workflow_rates.nc ${WFDIR}/wf_2e/tlim1.5win0.25/total-workflow_rates.nc --outfile ./tlim1.5win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p tlim2.0win0.25
rm ./tlim2.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25_GIS_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim2.0win0.25_globalsl_rates.nc ./tlim2.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim2.0win0.25/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim2.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim2.0win0.25/icesheets-pb1e-icesheets-tlim2.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/tlim2.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2e/tlim2.0win0.25/total-workflow_rates.nc --outfile ./tlim2.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p tlim3.0win0.25
rm ./tlim3.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25_GIS_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim3.0win0.25_globalsl_rates.nc ./tlim3.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim3.0win0.25/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim3.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim3.0win0.25/icesheets-pb1e-icesheets-tlim3.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/tlim3.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2e/tlim3.0win0.25/total-workflow_rates.nc --outfile ./tlim3.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p tlim4.0win0.25
rm ./tlim4.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25_GIS_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim4.0win0.25_globalsl_rates.nc ./tlim4.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim4.0win0.25/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim4.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim4.0win0.25/icesheets-pb1e-icesheets-tlim4.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/tlim4.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2e/tlim4.0win0.25/total-workflow_rates.nc --outfile ./tlim4.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10

# Remove the old files
mkdir -p tlim5.0win0.25
rm ./tlim5.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25_GIS_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim5.0win0.25_globalsl_rates.nc ./tlim5.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl_rates.nc ./tlim5.0win0.25/


# Generate the pbox component and total files
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25_AIS_globalsl_rates.nc ${PROJDIR}/icesheets-ipccar6-larmipicesheet-tlim5.0win0.25_TOT_globalsl_rates.nc --outfile ./tlim5.0win0.25/icesheets-pb1e-icesheets-tlim5.0win0.25_AIS_globalsl_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
srun python ${SCRIPTDIR}/generate_pbox_components_dev.py --infiles ${WFDIR}/wf_1e/tlim5.0win0.25/total-workflow_rates.nc ${WFDIR}/wf_2e/tlim5.0win0.25/total-workflow_rates.nc --outfile ./tlim5.0win0.25/total-workflow_rates.nc --pyear_start 2020 --pyear_end 2100 --pyear_step 10
