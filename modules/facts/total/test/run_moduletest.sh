#!/bin/bash
TESTSCRIPT_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $TESTSCRIPT_RELATIVE_DIR
TESTSCRIPT_DIR=`pwd`
FACTSSCRIPT_DIR=$TESTSCRIPT_DIR/../../../scripts
LOG_STDOUT=$TESTSCRIPT_DIR/test.out
LOG_STDERR=$TESTSCRIPT_DIR/test.err
PRESERVE_WORKDIR0=$PRESERVE_WORKDIR
SKIP_COLLECTION=1

source $FACTSSCRIPT_DIR/moduletest/moduletest.config.global
source moduletest.config

echo "Launching test script. Logging output to $LOG_STDOUT and $LOG_STDERR." 
source $FACTSSCRIPT_DIR/moduletest/moduletest.sh > $LOG_STDOUT 2> $LOG_STDERR

cp $WORKDIR/$GLOBALDIR/total-workflow.nc $OUTPUT_DIR/total-workflow-globalsl.nc
cp $WORKDIR/$LOCALDIR/total-workflow.nc $OUTPUT_DIR/total-workflow-localsl.nc

PRESERVE_WORKDIR=$PRESERVE_WORKDIR0
unset SKIP_COLLECTION

if [ -z "$PRESERVE_WORKDIR" ]
then
    echo ""
    echo "Cleaning up..."
    rm -fr $WORKDIR
else
    echo ""
    echo "Skipping clean up..."
fi


