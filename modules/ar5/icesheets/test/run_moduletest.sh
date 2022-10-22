#!/bin/bash
TESTSCRIPT_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $TESTSCRIPT_RELATIVE_DIR
TESTSCRIPT_DIR=`pwd`
FACTSSCRIPT_DIR=$TESTSCRIPT_DIR/../../../../scripts

source $FACTSSCRIPT_DIR/moduletest/moduletest.config.global
source moduletest.config
source $FACTSSCRIPT_DIR/moduletest/moduletest.sh

