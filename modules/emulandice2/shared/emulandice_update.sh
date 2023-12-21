#!/bin/bash

git submodule update --remote
Rscript -e "packrat::install('emulandice2')"
tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
