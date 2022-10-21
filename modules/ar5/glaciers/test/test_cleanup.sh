#!/bin/bash

mkdir $SCRIPT_DIR/output_global
mkdir $SCRIPT_DIR/output_local
mv *globalsl* $SCRIPT_DIR/output_global
mv *localsl* $SCRIPT_DIR/output_local

cd $SCRIPT_DIR
rm -fr $WORKDIR

