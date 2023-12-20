#!/bin/bash

git submodule update
Rscript emulandice_bundle_dependencies.R
tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
