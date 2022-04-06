#!/bin/bash

# Workflow and projection component directories --------------------------------
PROJDIR="/projects/kopp/ar6/global/full_sample_components"
WFDIR="/projects/kopp/ar6/global/full_sample_workflows/wf_3f"

# Go to the workflow directory
cd ${WFDIR}

# Remove the old files
mkdir -p ssp126
rm ./ssp126/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp126_GIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp126_AIS_globalsl.nc ./ssp126/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp126_globalsl.nc ./ssp126/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp126_globalsl.nc ./ssp126/


# Remove the old files
mkdir -p ssp245
rm ./ssp245/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp245_GIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp245_AIS_globalsl.nc ./ssp245/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp245_globalsl.nc ./ssp245/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp245_globalsl.nc ./ssp245/


# Remove the old files
mkdir -p ssp585
rm ./ssp585/*.nc

# Copy the new ones
cp ${PROJDIR}/glaciers-ar5-glaciersgmip2-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-FittedISMIP-icesheets-ssp585_GIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/icesheets-dp20-icesheet-ssp585_AIS_globalsl.nc ./ssp585/
cp ${PROJDIR}/oceandynamics-tlm-oceandynamics-ssp585_globalsl.nc ./ssp585/
cp ${PROJDIR}/landwaterstorage-ssp-landwaterstorage-ssp585_globalsl.nc ./ssp585/
