#!/bin/bash

git submodule update --remote
Rscript -e "packrat::install('emulandice2')"
