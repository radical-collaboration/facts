.. _chapter_quickstart:

Quick Start
===========

FACTS uses `RADICAL-EnsembleToolkit (EnTK) <https://radicalentk.readthedocs.io/en/stable/>`_ and `RADICAL-Pilot <https://radicalpilot.readthedocs.io/en/stable/>`_ to execute its modules. While the RADICAL tools are specifically designed for executing on a range of `supported <https://radicalpilot.readthedocs.io/en/stable/supported.html>`_ high performance computing (HPC) platforms, FACTS can also execute on a GNU/Linux workstation, virtual machine or container. Here we offer instructions for each deployment scenario but if you want to run FACTS at scale on an HPC platform, please `contact us <https://github.com/radical-collaboration/facts/issues/new>`_ and we will be happy to offer tailored support.

.. note:: Starting from version 1.4, the RADICAL tools will not require MongoDB anymore.

.. warning:: FACTS MUST be used within a dedicated Python virtual environment. You can use `venv`, `conda` or `virtualenv` to create one. If you try to install FACTS system-wide, it will fail.

Installing and Using FACTS on a GNU/Linux Workstation
-----------------------------------------------------

1. Clone the FACTS repository::

    git clone https://github.com/radical-collaboration/facts.git

2. Download modules-data.

   Archived versions are available on Zenodo at https://doi.org/10.5281/zenodo.7478191 and https://doi.org/10.5281/zenodo.7478447 (note, split between
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

4. Create and activate a Python virtual environment, and install FACTS's Python dependences in it. You can use `venv`, `conda` or `virtualenv` to create your Python virtual environment. See `these instructions <https://radicalpilot.readthedocs.io/en/stable/getting_started.html#Installation>`_ for further details. Using `venv`::

    python3 -m venv ve3
    . ve3/bin/activate
    pip install --upgrade setuptools pip wheel
    pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts
    pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq

5. Test your install by running the dummy experiment::

    python3 runFACTS.py experiments/dummy

6. Create a new experiment. For example::

    mkdir test
    cp -r experiments/coupling.ssp585/config.yml test

7. Run your experiment::

    python3 runFACTS.py test

Note that if you are running FACTS using localhost as a resource, all the input files for the experiment (which can be tens of GB) will get copied to ``~/radical.pilot.sandbox``. If you have space limits on your home directory, you may want to make this a symlink to a directory with fewer space limits prior to running FACTS.

If you wish to run the ``emulandice`` module set, additional steps are necessary, as this module set is a wrapper around separately developed R code (see https://github.com/tamsinedwards/emulandice/).

8. Ensure R and cmake are installed. On Ubuntu, these are provided by the r-base and cmake packages.

9. Build ``emulandice`` and a tar file of its associated R dependencies::

    modules/emulandice/emulandice_config.sh

Note that the data files for a FACTS experiment and transfered to the compute
resource with each experiment run. Thus, while it might in principle be possible
to run FACTS on your desktop and use a remote HPC resource, you probably don't
want to do this. At a minimum, you will want to have a fast, high-capacity
network connection to the resource.

Installing and Using FACTS on a GNU/Linux Virtual Machine or Container
----------------------------------------------------------------------

The RADICAL tools does not support MacOS or Windows. Therefore, to run on a Mac or Windows (the latter with WSL2), you need to run within a Linux virtual machine or container. 

.. warning:: These are not officially supported solutions; use them on your own recognizance.

To use a virtual machine on MacOS or Windows, you may want to investigate tools like `VirtualBox <https://www.virtualbox.org/>`_ or other commercial solutions. Once you create, run and log into a GNU/Linux VM, you can follow the instructions above to install and using FACTS.

Alternatively, you can use a Docker container. With Docker installed, you can launch an Ubuntu Focal environment::

    docker run --hostname=localhost --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch --runtime=runc -it ubuntu:focal.

This command assumes you have facts cloned into ``$HOME/facts`` and a writable scratch directory in ``$HOME/tmp``. Within the container , ``$HOME/facts`` will mount as ``/opt/facts`` and ``$HOME/tmp`` will mount as ``/scratch``.

Within this Ubuntu environment, the script `vm_factsenvsetup.sh <https://github.com/radical-collaboration/facts/blob/main/scripts/vm_factsenvsetup.sh>`_ will install and launch Mongo, install a suitable Python environment, install R and the dependencies of the ``emulandice`` module, and run the dummy experiment. This script may also be a helpful guide for installing FACTS in other clean environments.

Testing a module with a shell script
------------------------------------

In some cases, particularly during module development, it may be desirable to call
a FACTS module outside the EnTK framework. This can be done using an experimental
shell-script writing feature in runFACTS. Performance is not guaranteed, and
multi-module experiments are very likely not to work without customization, as
module coupling within FACTS is handled by the EnTK framework. 

1. Create an experiment (e.g., ``experiments/onemodule``) that invokes only the module of interest.

2. Create a shell scripts that executes the experiment by calling ``runFACTS`` with the ``--shellscript`` argument. For example::

    python3 runFACTS.py --shellscript experiments/onemodule > test.sh
    
3. Execute the shell script. For example::

    source test.sh
