#!/bin/sh

MODULETEST_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $MODULETEST_RELATIVE_DIR
MODULETEST_DIR=`pwd`
FACTSROOT_DIR=$MODULETEST_DIR/../..
TOTALOUTPUT_GLOBAL_DIR=$FACTSROOT_DIR/modules/total/test/data_global
TOTALOUTPUT_LOCAL_DIR=$FACTSROOT_DIR/modules/total/test/data_local

unset WORKDIR
$FACTSROOT_DIR/modules/kopp14/verticallandmotion/test/run_moduletest.sh

unset WORKDIR
$FACTSROOT_DIR/modules/ssp/landwaterstorage/test/run_moduletest.sh

unset WORKDIR
$FACTSROOT_DIR/modules/ipccar6/gmipemuglaciers/test/run_moduletest.sh

unset WORKDIR
$FACTSROOT_DIR/modules/ipccar6/ismipemuicesheet/test/run_moduletest.sh

unset WORKDIR
$FACTSROOT_DIR/modules/tlm/oceandynamics/test/run_moduletest.sh

if ! [[ -d $TOTALOUTPUT_GLOBAL_DIR ]]; then
    mkdir $TOTALOUTPUT_GLOBAL_DIR
else
    rm $TOTALOUTPUT_GLOBAL_DIR/*.nc
fi
cd $TOTALOUTPUT_GLOBAL_DIR

ln -s ../../../ssp/landwaterstorage/test/output/*globalsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*WAIS_globalsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*EAIS_globalsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*GIS_globalsl.nc
ln -s ../../../ipccar6/gmipemuglaciers/test/output/*globalsl.nc
ln -s ../../../tlm/oceandynamics/test/output/*globalsl.nc

if ! [[ -d $TOTALOUTPUT_LOCAL_DIR ]]; then
    mkdir $TOTALOUTPUT_LOCAL_DIR
else
    rm $TOTALOUTPUT_LOCAL_DIR/*.nc
fi
cd $TOTALOUTPUT_LOCAL_DIR

ln -s ../../../kopp14/verticallandmotion/test/output/*localsl.nc .
ln -s ../../../ssp/landwaterstorage/test/output/*localsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*WAIS_localsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*EAIS_localsl.nc
ln -s ../../../ipccar6/ismipemuicesheet/test/output/*GIS_localsl.nc
ln -s ../../../ipccar6/gmipemuglaciers/test/output/*localsl.nc
ln -s ../../../tlm/oceandynamics/test/output/*localsl.nc
