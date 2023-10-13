#!/bin/bash

#wget https://cran.r-project.org/src/contrib/Archive/dummies/dummies_1.5.6.tar.gz
Rscript emulandice_bundle_dependencies.R
tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
