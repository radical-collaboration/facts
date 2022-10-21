#!/bin/sh
SCRIPT_RELATIVE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd $SCRIPT_RELATIVE_DIR
SCRIPT_DIR=`pwd`

source $SCRIPT_DIR/test_init.sh
source $SCRIPT_DIR/test_run.sh
source $SCRIPT_DIR/test_cleanup.sh
