# Framework for Assessing Changes To Sea-level (FACTS) 

We are still catching up on things since our work with the Intergovernmental Panel on Climate Change's Sixth Assessment Report (IPCC AR6). This involves cleaning up lots of framework and module code. Your patience during this process is greatly appreciated!

See [Zenodo](https://doi.org/10.5281/zenodo.6419954) for the FACTS modules, data sets, and scripts used to produce the AR6 sea level projections.

See the [IPCC-AR6-Sea-Level-Projections repo](https://github.com/rutgers-ESSP/IPCC-AR6-Sea-Level-Projections) for a guide for accessing the IPCC AR6 sea level projections data. 

### Module Tests

Almost all modules have test scripts that allow them to be run outside the FACTS/EnTK framework. The configuration of the module test scripts are specified in a moduletest.config file. See, for example, [modules/ar5/icesheets/test/moduletest.config](modules/ar5/icesheets/test/moduletest.config). There may also be global settings (e.g., the scratch directory you want used) that need to be set in [scripts/moduletest/moduletest.config.global](scripts/moduletest/moduletest.config.global).

Since, in running modules outside the FACTS/EnTK framework, you will not have the benefits of EnTK's environment management, you will need to make sure all the packages needed to support the modules are installed in their environment. This will differ between packages (e.g., [emulandice](modules/emulandice) is a FACTS wrapper around independently developed R code, and running it requires all the R packages required by that code), but a good working environment for most purposes can be set up with conda as follows:

```
conda create -n facts python=3.7 -c conda-forge -y
conda activate facts
conda install radical.entk -c conda-forge
conda install numpy scipy netCDF4 pyyaml matplotlib h5py yq
```

In addition, you will need to have all the associated module data tgz files in your [module-data](module-data) directory.