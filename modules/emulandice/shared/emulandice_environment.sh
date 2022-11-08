#! /bin/bash

# this script loads the modules
# required for emulandice to run on Amarel
#
# You will have to modify it for your system.

module use /projects/community/modulefiles
module load intel/19.0.3
module load R/3.6.3-gc563
module load gcc/5.4

if [ "$1" == "--Rscript" ]; then
    Rscript -e "source('packratc/init.R')" -e "packrat::install_local('emulandice_1.1.0.tar.gz')"
fi
