.. _chapter_quickstart:

Quick Start
===========

Installing and Using FACTS
--------------------------

1. Clone the FACTS repository::

    git clone https://github.com/radical-collaboration/facts.git

2. Download modules-data.

   Archived versions are available on Zenodo at doi:10.5281/zenodo.7478191 and doi:10.5281/zenodo.7478447 (note, split between
   two Zenodo entries because of size limitations), while a development version is currently synced at 
   https://rutgers.box.com/s/6vjio67b533lx5vzgyt5e5lw5jb4ftts. (If you have multiple users of FACTS, you might want to put
   these ~60 GB of files in a common location and soft-link to each user's directory.)

3. Set up and launch MongoDB server. Options include:

  - `Install MongoDB Community Edition on Linux <https://www.mongodb.com/docs/manual/administration/install-on-linux/>`_

  - Run MongoDB from a container. On a system with Singularity installed, this looks something like::

      mkdir mongo
      singularity build --sandbox mongo/ docker://mongo
      singularity run -w mongo &.

  - Set up your resource file to use the MongoDB server run by RADICAL. Ask for MongoDB parameters by writing to the FACTS
    team via email or by opening an issue in this repository.

4. Create and activate a Python virtual environment, and install FACTS's Python
   dependences in it. Using `venv`::

    python3 -m venv ve3
    . ve3/bin/activate
    pip install --upgrade setuptools pip wheel
    pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts
    pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq.

5. Test your install by running the dummy experiment::

    python3 runFACTS.py experiments/dummy

6. Create a new experiment. For example::

    mkdir test
    cp -r experiments/coupling.ssp585/config.yml experiments/coupling.ssp585/locations.lst test

7. Run your experiment::

    python3 runFACTS.py test

Note that if you are running FACTS using localhost as a resource, all the input
files for the experiment (which can be tens of GB) will get copied to
```~/radical.pilot.sandbox```. If you have space limits on your home directory,
you may therefore want to make this a symlink to a directory with fewer space
limits prior to running.

Using FACTS on a Mac
--------------------

The RADICAL stack does not support MacOS. Therefore, to run on a Mac, you need
to run within a Linux environment. One way to do this is with a Docker
container. Note this is not an officially supported solution, and you must do so
on your own recognizance.

With Docker installed, you can launch an Ubuntu environment::

    docker run --hostname=localhost --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch --runtime=runc -it ubuntu:focal.

This command assumes you have facts installed in ```$HOME/facts``` and a
writable scratch directory in ```$HOME/tmp```.

Within this Ubuntu environment, the script mac_docker_factsenvsetup.sh will
install Mongo and a suitable Python environment, and run the dummy experiment.

This solution may also work on Windows, but has not been tested.

Note that the data files for a FACTS experiment and transfered to the compute
resource with each experiment run. Thus, while it might in principle be possible
to run FACTS on your desktop and use a remote HPC resource, you probably don't
want to do this. At a minimum, you will want to have a fast, high-capacity
network connection to the resource.

Testing a module with a shell script
------------------------------------

In some cases, it may be desirable to call a FACTS module outside the EnTK framework.
This can be done using an experimental shell-script writing feature in runFACTS.
Performance is not guaranteed, and multi-module experiments are very likely not to
work without customization. 

1. Create an experiment (e.g., ``experiments/onemodule``) that invokes only the module of interest.

2. Create a shell scripts that executes the experiment by calling ``runFACTS`` with the ``--shellscript`` argument. For example::

    python3 runFACTS.py --shellscript experiments/onemodule > test.sh
    
3. Execute the shell script. For example::

    source test.sh
