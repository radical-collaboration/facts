#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/full_sample_components"
WFDIR="/projects/kopp/ar6/global/full_sample_workflows/wf_4"

# Go to the workflow directory
cd ${WFDIR}

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_AIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp126_GIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl.nc ./ssp126/


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_AIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp245_GIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./ssp245/


# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_AIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-ssp585_GIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl.nc ./ssp585/


# Remove the old files
mkdir -p tlim2.0win0.25
rm ./tlim2.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim2.0win0.25_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-tlim2.0win0.25_AIS_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-tlim2.0win0.25_GIS_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim2.0win0.25_globalsl.nc ./tlim2.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim2.0win0.25/


# Remove the old files
mkdir -p tlim5.0win0.25
rm ./tlim5.0win0.25/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-tlim5.0win0.25_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-tlim5.0win0.25_AIS_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/icesheets-ipccar6-bambericesheet-tlim5.0win0.25_GIS_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-tlim5.0win0.25_globalsl.nc ./tlim5.0win0.25/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./tlim5.0win0.25/
