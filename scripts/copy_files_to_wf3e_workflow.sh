#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/regional/full_sample_components"
WFDIR="/projects/kopp/ar6/regional/full_sample_workflows/wf_3e"

# Go to the workflow directory
cd ${WFDIR}

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp126_GIS_localsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp126_AIS_localsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_localsl.nc ./ssp126/
cp ${PROJDIR}/verticallandmotion-kopp14-verticallandmotion_localsl.nc ./ssp126/

# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp245_GIS_localsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp245_AIS_localsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_localsl.nc ./ssp245/
cp ${PROJDIR}/verticallandmotion-kopp14-verticallandmotion_localsl.nc ./ssp245/

# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ipccar6-gmipemuglaciers-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-ismipemuicesheet-ssp585_GIS_localsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp585_AIS_localsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_localsl.nc ./ssp585/
cp ${PROJDIR}/verticallandmotion-kopp14-verticallandmotion_localsl.nc ./ssp585/
