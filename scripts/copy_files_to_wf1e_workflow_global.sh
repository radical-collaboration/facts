#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/full_sample_components"
WFDIR="/projects/kopp/ar6/global/full_sample_workflows/wf_1e"

# Go to the workflow directory
cd ${WFDIR}

# Remove the old files
mkdir -p ssp119
rm ./ssp119/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp119_globalsl.nc ./ssp119/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp119_GIS_globalsl.nc ./ssp119/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp119_AIS_globalsl.nc ./ssp119/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp119_globalsl.nc ./ssp119/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp119_globalsl.nc ./ssp119/


# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_GIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_AIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl.nc ./ssp126/


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_GIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_AIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./ssp245/


# Remove the old files
mkdir -p ssp370
rm ./ssp370/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp370_globalsl.nc ./ssp370/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp370_GIS_globalsl.nc ./ssp370/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp370_AIS_globalsl.nc ./ssp370/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp370_globalsl.nc ./ssp370/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp370_globalsl.nc ./ssp370/


# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_GIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_AIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl.nc ./ssp585/


# Remove the old files
mkdir -p tlim1.5win0.25
rm ./tlim1.5win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim1.5win0.25_globalsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25_GIS_globalsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim1.5win0.25_AIS_globalsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim1.5win0.25_globalsl.nc ./tlim1.5win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim1.5win0.25/


# Remove the old files
mkdir -p tlim2.0win0.25
rm ./tlim2.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim2.0win0.25_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25_GIS_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim2.0win0.25_AIS_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim2.0win0.25_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim2.0win0.25/


# Remove the old files
mkdir -p tlim3.0win0.25
rm ./tlim3.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim3.0win0.25_globalsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25_GIS_globalsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim3.0win0.25_AIS_globalsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim3.0win0.25_globalsl.nc ./tlim3.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim3.0win0.25/


# Remove the old files
mkdir -p tlim4.0win0.25
rm ./tlim4.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim4.0win0.25_globalsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25_GIS_globalsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim4.0win0.25_AIS_globalsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim4.0win0.25_globalsl.nc ./tlim4.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim4.0win0.25/


# Remove the old files
mkdir -p tlim5.0win0.25
rm ./tlim5.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-tlim5.0win0.25_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25_GIS_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-tlim5.0win0.25_AIS_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim5.0win0.25_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim5.0win0.25/
