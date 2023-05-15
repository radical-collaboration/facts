# emulandice

The emulandice family of modules wrap around the Edwards et al. (2021) Gaussian process emulators of the ISMIP6 and GlacierMIP2 models.

Before running any emulandice module you should run
```
emulandice_config.sh
```
to create two files needed by the modules within the emulandice package, emulandice_1.1.0.tar.gz and emulandice_bundled_dependencies.tgz. 

You may need to customize shared/emulandice_environment.sh and shared/emulandice_bundle_dependencies.R based on your local environment.

If you have not run R before and do not have a user-writeable directory into which R packages can be installed, you may want to manually open R and run install.packages("packrat").

