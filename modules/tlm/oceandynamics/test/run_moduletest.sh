#!/bin/bash
TESTSCRIPT_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $TESTSCRIPT_RELATIVE_DIR
TESTSCRIPT_DIR=`pwd`
FACTSSCRIPT_DIR=$TESTSCRIPT_DIR/../../../../scripts
LOG_STDOUT=$TESTSCRIPT_DIR/test.out
LOG_STDERR=$TESTSCRIPT_DIR/test.err


source $FACTSSCRIPT_DIR/moduletest/moduletest.config.global
source moduletest.config

echo "Launching test script. Logging output to $LOG_STDOUT and $LOG_STDERR." 
source $FACTSSCRIPT_DIR/moduletest/moduletest.sh > $LOG_STDOUT 2> $LOG_STDERR


