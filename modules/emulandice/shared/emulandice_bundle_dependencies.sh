#!/bin/bash

Rscript emulandice_bundle_dependencies.R
tar cvzf emulandice_bundled_dependencies.tgz .Rprofile packrat/
