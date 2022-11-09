#!/bin/bash

# This script is intended to be used to support
# the development of scripts that call
# modules for testing purposes outside the FACTS EnTK framework.

echo "Initiating test in $WORKDIR..."

if [ -z "$SKIP_WORKDIR_SETUP" ]; then  
    mkdir $WORKDIR
    cd $WORKDIR
    cp -L -r $TESTSCRIPT_DIR/../* .
else
    cd $WORKDIR
fi

if [ -z "$SKIP_EXTRACT" ]
then
    echo ""
    echo "Extracting data files..."
    for i in data/*
    do
        ln -s $i .
    done
    for i in data/*.tgz
    do
	tar xzf $i  2>&1 | grep -v 'Ignoring'
    done
else
    echo ""
    echo "Skipping data extraction..."
fi

echo ""
echo "Executing workflow..."

j=0
for i in "${PREEXEC[@]}"
do
    j=$(( $j+1 ))
    EXECCMD="${PREEXEC[j]}"
    echo $EXECCMD
    $EXECCMD
done

j=0
for i in "${STAGES[@]}"
do
    j=$(( $j+1 ))
    EXECCMD="python ${STAGE_SCRIPT[j]} ${STAGEOPTIONS[j]}"
    echo $EXECCMD
    $EXECCMD
done


echo ""
echo "Collecting output files..."
if ! [[ -d $OUTPUT_DIR ]]
then
    mkdir $OUTPUT_DIR
fi

if [ -z "$SKIP_COLLECTION" ]
then
    ls ${PIPELINE_ID}*.nc
    mv ${PIPELINE_ID}*.nc $OUTPUT_DIR
fi

cd $TESTSCRIPT_DIR

if [ -z "$PRESERVE_WORKDIR" ]
then
    echo ""
    echo "Cleaning up..."
    rm -fr $WORKDIR
else
    echo ""
    echo "Skipping clean up..."
fi
