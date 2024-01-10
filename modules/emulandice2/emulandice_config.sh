#!/bin/bash

ROOTDIR=`dirname $0`
MODE=$1


cd $ROOTDIR/shared
if [ -f emulandice2/emulandice2.Rproj ]; then
    echo "emulandice2 R project already cloned"
else
    echo "Initializing emulandice2 R project"
    git submodule update --init
fi

case $MODE in 
    "clean")
        echo Cleaning up old files...
        rm -v -f emulandice_*gz* emulandice_bundled_dependencies.tgz .Rprofile
        rm -v -fr packrat/
        ;;
    *)
        source emulandice_environment.sh

        FILE="emulandice_bundled_dependencies.tgz"
        if [ -f $FILE ]; then
            echo "$FILE already exists."
        else
            echo "Building $FILE..."
            source emulandice_bundle_dependencies.sh
        fi
        ;;
esac

echo
echo If you have not already, please download emulandice2\*
echo from https://rutgers.box.com/v/facts-module-data
echo and place in the modules-data directory.
