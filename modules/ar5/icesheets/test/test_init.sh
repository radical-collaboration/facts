#!/bin/bash

cd $SCRIPT_RELATIVE_DIR
SCRIPT_DIR=`pwd`

WORKDIR=/scratch/`whoami`/test.`date +%s`
mkdir $WORKDIR
cd $WORKDIR

pwd

cp -L -r $SCRIPT_DIR/../* .
LFILE=test/location_list.lst

for i in data/*
do
    tar xzvf $i  2>&1 | grep -v 'Ignoring'
done
