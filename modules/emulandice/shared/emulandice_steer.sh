#!/bin/bash
if [[ -d results ]]
    mkdir results
fi

Rscript -e "source('packrat/init.R')"

# have to do something special for dummies since it's been removed
Rscript -e 'packrat::set_opts(local.repos = "packrat/src")' -e 'packrat::install_local("dummies/dummies_1.5.6.tar.gz")'

# restore other packages
Rscript -e "source('packrat/init.R')" -e 'packrat::restore()'

# install emulandice package
Rscript -e 'packrat::set_opts(local.repos = "packrat/src/emulandice")' -e 'packrat::install_local("emulandice_1.1.0.tar.gz")'

# run emulandice
emulandice_dataset = $1
N_FACTS=$2
Rscript -e 'library(emulandice)' -e "main('decades',dataset=$emulandice_dataset,N_FACTS=$N_FACTS)"
