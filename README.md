# Framework for Assessing Changes To Sea-level (FACTS)

The Framework for Assessing Changes To Sea-level (FACTS) is an open-source modular, scalable framework for global mean, regional, and extreme sea level projection that is designed to support the characterization of ambiguity in sea-level projections. It is designed so users can easily explore deep uncertainty by investigating the implications on GMSL, RSL, and ESL of different choices for different processes. Its modularity allows components to be represented by either simple or complex model. 

## Installing and Using FACTS

1. Clone the FACTS repository:

  ```
  git clone git@github.com:radical-collaboration/facts.git
  ```

2. Download modules-data. 
 
   Archived versions are available on Zenodo at [doi:10.5281/zenodo.7478191](https://doi.org/10.5281/zenodo.7478191) and [doi:10.5281/zenodo.7478447](https://doi.org/10.5281/zenodo.7478447) (note, split between two Zenodo entries because of size limitations), while a development version is currently synced at https://rutgers.box.com/s/6vjio67b533lx5vzgyt5e5lw5jb4ftts. (If you have multiple users of FACTS, you might want to put these ~60 GB of files in a common location and soft-link to each user's directory.) 

3. Set up and launch MongoDB server. Options include:
 
    * Install MongoDB server. Follow the [official documentation](https://www.mongodb.com/docs/manual/administration/install-on-linux/).

    * Run MongoDB from a container. On a system with Singularity installed, this looks something like:
 
      ```
      mkdir mongo
      singularity build --sandbox mongo/ docker://mongo
      singularity run -w mongo &
      ```

    * Set up your resource file to use the MongoDB server run by RADICAL. Ask for MongoDB parameters by writing to the FACTS team via email or by opening an issue in this repository.

4. Create and activate a Python virtual environment, and install FACTS's Python dependences in it:

   - Using `venv`:

     ```
     python3 -m venv ve3
     . ve3/bin/activate
     pip install --upgrade setuptools pip wheel
     pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts
     pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq
     ```

5. Execute FACTS.

   ```
   mkdir test
   cp -r experiments/coupling.ssp585/config.yml experiments/coupling.ssp585/locations.lst test
   python3 runFACTS.py test
   ```

Note that if you are running FACTS using localhost as a resource, all the input files for the experiment (which can be tens of GB) will get copied to ```~/radical.pilot.sandbox```. If you have space limits on your home directory, you may therefore want to make this a symlink to a directory with fewer space limits prior to running.

## Using FACTS on a Mac

The RADICAL stack does not support MacOS. Therefore, to run on a Mac, you need to run within a Linux environment. One way to do this is with a Docker container. Note this is not an officially supported solution, and you must do so on your own recognizance.

With Docker installed, you can launch an Ubuntu environment:

  ```
  docker run --hostname=localhost --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch --runtime=runc -it ubuntu:focal
  ```
(This command assumes you have facts installed in ```$HOME/facts``` and a writable scratch directory in ```$HOME/tmp```.)

Within this Ubuntu environment, the script mac_docker_factsenvsetup.sh will install Mongo and a suitable Python environment, and run the dummy experiment.

This solution may also work on Windows, but has not been tested.

Note that the data files for a FACTS experiment and transfered to the compute resource with each experiment run. Thus, while it might in principle be possible to run FACTS on your desktop and use a remote HPC resource, you probably don't want to do this. At a minimum, you will want to have a fast, high-capacity network connection to the resource.

## Resources

See [doi:10.5281/zenodo.6419954](https://doi.org/10.5281/zenodo.6419954) for the FACTS modules, data sets, and scripts used to produce the AR6 sea level projections.

See the [IPCC-AR6-Sea-Level-Projections repo](https://github.com/rutgers-ESSP/IPCC-AR6-Sea-Level-Projections) for a guide for accessing the IPCC AR6 sea level projections data.

