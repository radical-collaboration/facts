#!/bin/bash

ROOTDIR=`dirname $0`
MODE=$1

case $MODE in 
    "clean")
	cd $ROOTDIR/shared
        echo Cleaning up old files...
        rm -v -f emulandice_*gz* emulandice_bundled_dependencies.tgz .Rprofile
        rm -v -fr packrat/
        ;;
    "update")
	if [ -f $ROOTDIR/shared/emulandice2/emulandice2.Rproj ]; then
		echo "Updating emulandice2 R project"
		git submodule update --remote
	else
    		echo "Initializing emulandice2 R project"
    		git submodule update --init
	fi

	cd $ROOTDIR/shared
	source emulandice_environment.sh

	Rscript -e "packrat::install('emulandice2')"
	tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
	;;
    "help")
        echo $0 clean: Remove bundled dependencies.
        echo $0 build: Build dependencies for emulandice2.
        echo $0 update: Update emulandice2 from repo and rebuild.
	;;
    *)
	if [ -f $ROOTDIR/shared/emulandice2/emulandice2.Rproj ]; then
		echo "emulandice2 R project already cloned"
	else
    		echo "Initializing emulandice2 R project"
    		git submodule update --init
	fi

	cd $ROOTDIR/shared
        source emulandice_environment.sh

        FILE="emulandice_bundled_dependencies.tgz"
        if [ -f $FILE ]; then
            echo "$FILE already exists."
        else
            echo "Building $FILE..."
	    Rscript emulandice_bundle_dependencies.R
	    tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
	fi

	echo
	echo If you have not already, please download emulandice2\*
	echo from https://rutgers.box.com/v/facts-module-data
	echo and place in the modules-data directory.
        ;;
esac

