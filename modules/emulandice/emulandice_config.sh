#!/bin/bash

ROOTDIR=`dirname $0`
MODE=$1

cd $ROOTDIR/shared
case $MODE in 
    "clean")
        echo Cleaning up old files...
        rm -v *.gz* *.tgz
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
