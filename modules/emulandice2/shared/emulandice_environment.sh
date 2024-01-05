#! /bin/bash

# this script loads the modules
# required for emulandice to run on Amarel
#
# You may have to modify it for your system.

hs=`hostname`

if [ "${hs: -18}" = "amarel.rutgers.edu" ]; then
    module use /projects/community/modulefiles

    module purge
    module load intel/19.0.3
    module load R/3.6.3-gc563
    module load gcc/7.3.0-gc563
    module load cmake/3.24.3-sw1088

    module list
fi

