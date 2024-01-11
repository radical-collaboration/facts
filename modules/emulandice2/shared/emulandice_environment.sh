#! /bin/bash

# this script loads the modules
# required for emulandice to run on Amarel
#
# You may have to modify it for your system.

hs=`hostname`

if [ "${hs: -18}" = "amarel.rutgers.edu" ]; then
    module use /projects/community/modulefiles

    module purge
    module load netcdf/4.9.0-sw1088
    module load R/4.1.0-gc563
    module load gcc/10.2.0-bz186
    module load cmake/3.24.3-sw1088

    module list
fi

