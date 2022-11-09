Before running the emulandice module, you should
```
source emulandice_environment.sh
source emulandice_build.sh
source emulandice_bundle_dependencies.sh
```
to create two files needed by the other modules within the emulandice package, emulandice_1.1.0.tar.gz and emulandice_bundled_dependencies.tgz. 

You will likely need to customize emulandice_environment.sh and emulandice_bundle_dependencies.R based on your local environment.

If you have not run R before and do not have a user-writeable directory into which R packages can be installed, you may want to manually open R and run install.packages("packrat").
