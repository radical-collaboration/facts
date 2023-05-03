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

        FILE="emulandice_1.1.0.tar.gz"
        if [ -f $FILE ]; then
            echo "$FILE alreadu exists."
        else
            echo "Building $FILE..."
            source emulandice_build.sh
        fi

        FILE="emulandice_bundled_dependencies.tgz"
        if [ -f $FILE ]; then
            echo "$FILE alreadu exists."
        else
            echo "Building $FILE..."
            source emulandice_bundle_dependencies.sh
        fi
        ;;
esac