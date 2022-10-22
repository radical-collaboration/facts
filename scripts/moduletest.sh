#!/bin/bash

# This script is intended to be used to support
# the development of scripts that call
# modules for testing purposes outside the FACTS EnTK framework.

echo "Initiating test in $WORKDIR..."

mkdir $WORKDIR
cd $WORKDIR

pwd

cp -L -r $TESTSCRIPT_DIR/../* .

echo "Extracting data files..."
for i in data/*
do
    tar xzvf $i  2>&1 | grep -v 'Ignoring'
done

echo "Executing workflow..."

STAGE_SCRIPTS=()
j=0
for i in $STAGES
do
    j=$(( j+1 ))
    STAGE_SCRIPT=(`yq .\"$i\".task1.script ../pipeline.yml`)
    python $STAGE_SCRIPT ${STAGEOPTIONS[j]}
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
rm -fr $WORKDIR
