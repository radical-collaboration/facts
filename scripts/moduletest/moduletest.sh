#!/bin/bash

# This script is intended to be used to support
# the development of scripts that call
# modules for testing purposes outside the FACTS EnTK framework.

echo "Initiating test in $WORKDIR..."

mkdir $WORKDIR
cd $WORKDIR

cp -L -r $TESTSCRIPT_DIR/../* .

if [[ ! -v SKIP_EXTRACTDATA ]]; then
    echo "Extracting data files..."
    for i in data/*
    do
	tar xzf $i  2>&1 | grep -v 'Ignoring'
    done
fi

echo "Executing workflow..."

j=0
for i in "${STAGES[@]}"
do
    j=$(( $j+1 ))
    EXECCMD="python ${STAGE_SCRIPT[j]} ${STAGEOPTIONS[j]}"
    echo $EXECCMD
    $EXECCMD
done


echo "Collecting output files..."
if ! [[ -d $OUTPUTGLOBAL_DIR ]]
then
    mkdir $OUTPUTGLOBAL_DIR
fi

if ! [[ -d $OUTPUTLOCAL_DIR ]]
then
    mkdir $OUTPUTLOCAL_DIR
fi

mv *globalsl* $OUTPUTGLOBAL_DIR
mv *localsl* $OUTPUTLOCAL_DIR

cd $TESTSCRIPT_DIR

if [[ ! -v PRESERVE_WORKDIR ]]; then
    rm -fr $WORKDIR
fi
