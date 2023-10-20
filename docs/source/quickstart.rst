.. _chapter_quickstart:

Quick Start
===========

FACTS uses `RADICAL-EnsembleToolkit (EnTK) <https://radicalentk.readthedocs.io/en/stable/>`_ and `RADICAL-Pilot <https://radicalpilot.readthedocs.io/en/stable/>`_ to execute its modules. While the RADICAL tools are specifically designed for executing on a range of `supported <https://radicalpilot.readthedocs.io/en/stable/supported.html>`_ high performance computing (HPC) platforms, FACTS can also execute on a GNU/Linux workstation, virtual machine or container. Here we offer instructions for each deployment scenario but if you want to run FACTS at scale on an HPC platform, please `contact us <https://github.com/radical-collaboration/facts/issues/new>`_ and we will be happy to offer tailored support.

.. note:: Starting from version 1.40, the RADICAL tools will not require MongoDB anymore.

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

3. Create and activate a Python virtual environment, and install FACTS's Python dependences in it. You can use `venv`, `conda` or `virtualenv` to create your Python virtual environment. See `these instructions <https://radicalpilot.readthedocs.io/en/stable/getting_started.html#Installation>`_ for further details. Using `venv`::

    python3 -m venv ve3
    . ve3/bin/activate
    pip install --upgrade setuptools pip wheel radical.entk pyyaml

4. Test your install by running the dummy experiment::

    python3 runFACTS.py experiments/dummy

5. If you wish to run the ``emulandice`` module set, additional steps are necessary, as this module set is a wrapper around separately developed R code (see https://github.com/tamsinedwards/emulandice/). First, ensure R and cmake are installed. On Ubuntu, these are provided by the r-base and cmake packages. Then build ``emulandice`` and a tar file of its associated R dependencies::

    modules/emulandice/emulandice_config.sh

6. Create a new experiment. For example::

    mkdir test
    cp -r experiments/coupling.ssp585/config.yml test

7. Run your experiment::

    python3 runFACTS.py test

Note that if you are running FACTS using localhost as a resource, all the input files for the experiment (which can be tens of GB) will get copied to ``~/radical.pilot.sandbox``. If you have space limits on your home directory, you may want to make this a symlink to a directory with fewer space limits prior to running FACTS.

Note that the data files for a FACTS experiment and transfered to the compute
resource with each experiment run. Thus, while it might in principle be possible
to run FACTS on your desktop and use a remote HPC resource, you probably don't
want to do this. At a minimum, you will want to have a fast, high-capacity
network connection to the resource.

Installing and Using FACTS on a GNU/Linux Container
----------------------------------------------------------------------

The RADICAL toolkit does not support MacOS or Windows. Therefore, to run on a Mac or Windows (the latter with WSL2), you need to run within a Linux virtual machine or container. 

We have provided a Docker container in the ``docker/`` directory. This container provides the Linux,
Python, R, and RADICAL toolkit environment needed for FACTS to run.
FACTS itself does not reside within the container because of needs related to
storage space for module data, persistence of changes, and writability. The instructions below
assume FACTS resides outside the container in ``$HOME/facts`` and mount it within the container as
``/opt/facts``. At the moment, the docker environment appears to work fairly reliably when
using localhost as the resource, but not when using remote resources. 

To install FACTS through Docker please follow the steps below:

1. Clone the FACTS repository::

    git clone https://github.com/radical-collaboration/facts.git

2. Download modules-data.

   Archived versions are available on Zenodo at https://doi.org/10.5281/zenodo.7478191 and https://doi.org/10.5281/zenodo.7478447 (note, split between
   two Zenodo entries because of size limitations), while a development version is currently synced at 
   https://rutgers.box.com/s/6vjio67b533lx5vzgyt5e5lw5jb4ftts.

3. cd into the ``docker`` directory

4. Build the docker container::

    sh develop.sh

5. Start a container from the ``facts`` image, assuming that the FACTS repository was cloned in ``$HOME/facts``::

    docker run -it --volume=$HOME/facts:/opt/facts -w /opt/facts facts

7. Confirm that FACTS work within the container::

    python3 runFACTS.py experiments/dummy

8. If you wish to use ``emulandice``, build ``emulandice`` and a tar file of its associated R dependencies::

    modules/emulandice/emulandice_config.sh

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
