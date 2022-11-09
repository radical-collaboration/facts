# Framework for Assessing Changes To Sea-level (FACTS)

The Framework for Assessing Changes To Sea-level (FACTS) is an open-source modular, scalable, and extensive framework for global mean, regional, and extreme sea level projection that is designed to support the characterization of ambiguity in sea-level projections. It is designed so users can easily explore deep uncertainty by investigating the implications on GMSL, RSL, and ESL of different choices for different processes. Its modularity allows components to be represented by either simple or complex model. Because it is built upon the Radical-PILOT computing stack, different modules can be dispatched for execution on resources appropriate to their computational complexity.

## Resources

See [Zenodo](https://doi.org/10.5281/zenodo.6419954) for the FACTS modules, data sets, and scripts used to produce the AR6 sea level projections.

See the [IPCC-AR6-Sea-Level-Projections repo](https://github.com/rutgers-ESSP/IPCC-AR6-Sea-Level-Projections) for a guide for accessing the IPCC AR6 sea level projections data.

## Using FACTS on a local Linux machine

1. Cloning FACTS repository on the local Linux machine:

  ```
  git clone git@github.com:radical-collaboration/facts.git
  ```

2. Download modules-data from https://rutgers.box.com/s/6vjio67b533lx5vzgyt5e5lw5jb4ftts. (If you have multiple users of FACTS, you might want to put these >50 GB of file in a common location and soft-link to each user's directory.) 

3. Installing MongoDB and RabbitMQ servers. Follow the official documentation:

  - [Install MongoDB Community Edition on Linux](https://www.mongodb.com/docs/manual/administration/install-on-linux/)
  - [Downloading and Installing RabbitMQ](https://www.rabbitmq.com/download.html)

4. Creating and activating a Python virtual environment, and installing FACTS's Python dependences in it:

  - Using `venv`:

    ```
    python3 -m venv ve3
    . ve3/bin/activate
    pip install --upgrade setuptools pip wheel
    pip install radical.entk
    pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq
    ```

  - Using `Conda`:

    ```
    conda create -n facts
    conda activate facts
    conda install -c conda-forge radical.entk
    conda install numpy scipy netCDF4 pyyaml matplotlib h5py yq
    ```

5. executing FACTS.

  ```
  mkdir test
  cp -r experiments/coupling.ssp585/* test
  python3 runFACTS.py test
  ```

## Using FACTS on Rutgers' Amarel

1. Login into Amarel, see [Amarel User Guide](https://oarc.rutgers.edu/resources/amarel/).

2. Cloning FACTS repository on the local Linux machine:

  ```
  git clone git@github.com:radical-collaboration/facts.git
  ```

3. Download modules-data from https://rutgers.box.com/s/6vjio67b533lx5vzgyt5e5lw5jb4ftts. If you have multiple users of FACTS, you might want to put these >50 GB of file in a common location and soft-link to each user's directory, for example:

```
cd ~/facts/modules-data
ln -s /projects/kopp/facts-dev/modules-data/*.tgz . 
```

4. Configuring MongoDB, RabbitMQ and Amarel resource parameters:

  - Ask for MongoDB and RabbitMQ parameters to the FACTS team via email or by opening an issue in this repository.
  - Edit FACTS configuration file:

    ```
    vim experiments/temp_exp/config.yml
    -> replace 'localhost' with 'amarel'
    -> enter the parameters for RabbitMQ/MongoDB
    ```

5. Creating and activating a Python virtual environment, and installing FACTS's Python dependences in it:

  - Using `venv`:

    ```
    python3 -m venv ve3
    . ve3/bin/activate
    pip install --upgrade setuptools pip wheel
    pip install radical.entk
    pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq
    ```

  - Using `Conda`:

    ```
    conda create -n facts
    conda activate facts
    conda install -c conda-forge radical.entk
    conda install numpy scipy netCDF4 pyyaml matplotlib h5py yq
    ```

6. executing FACTS.

  ```
  mkdir test
  cp -r experiments/coupling.ssp585/* test
  python3 runFACTS.py test
  ```


### Module Tests

Almost all modules have test scripts that allow them to be run outside the FACTS/EnTK framework. These should be invoked via the test/run_moduletest.sh script. The configuration of the module test scripts are specified in a moduletest.config file. See, for example, [modules/ar5/icesheets/test/moduletest.config](modules/ar5/icesheets/test/moduletest.config). There may also be global settings (e.g., the scratch directory you want used) that need to be set in [scripts/moduletest/moduletest.config.global](scripts/moduletest/moduletest.config.global).

Since, in running modules outside the FACTS/EnTK framework, you will not have the benefits of EnTK's environment management, you will need to make sure all the packages needed to support the modules are installed in their environment. This will differ between packages (e.g., [emulandice](modules/emulandice) is a FACTS wrapper around independently developed R code, and running it requires all the R packages required by that code), but a good working environment for most purposes can be set up with conda as follows:

```
conda create -n facts python=3.7 -c conda-forge -y
conda activate facts
conda install radical.entk -c conda-forge
conda install numpy scipy netCDF4 pyyaml matplotlib h5py yq
```

In addition, you will need to have all the associated module data tgz files in your [modules-data](modules-data) directory. Note that some of these files are quite large (the total exceeds 50 GB), so if you have multiple users on a system employing FACTS, best practice would be to have a shared directory in which all these large files live and then use sym-links to link them to each user's modules-data directory, e.g.:

```
cd ~/facts/modules-data
ln -s /projects/shared/facts/modules-data/*.tgz .
```
