#!/bin/bash

# This script will create an environment on a virtual machine that can run FACTS.
#
# At least on an arm64 Mac, it will work in docker starting with the command:
#
#    docker run --hostname=localhost --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch --runtime=runc -it ubuntu:focal

# update these two lines for your system
FACTSROOT=/opt/facts # the location of the FACTS main directory
RADICALPILOTSANDBOX=/scratch/radical.pilot.sandbox # where you want the radical pilot sandbox to live -- if empty, will use home dir

apt-get update -y

# Install MonogDB following directions at
# https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/
# (you will need to update the line amending the apt sources list if running on a distribution other than Ubuntu Focal)

apt-get install -y sudo gnupg curl apt-transport-https
curl -fsSL https://pgp.mongodb.com/server-6.0.asc | \
   sudo gpg -o /usr/share/keyrings/mongodb-server-6.0.gpg \
   --dearmor
echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/6.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-6.0.list
apt-get update -y
apt-get install -y mongodb-org

# Start Mongod
/usr/bin/mongod --config /etc/mongod.conf --fork

## Install python
apt-get install -y python3-pip python3-venv git libnetcdf-dev python3-netcdf4
apt-get install -y ssh iputils-ping slurm

# Set up Python environment
python3 -m venv ve3
. ve3/bin/activate
pip install --upgrade setuptools pip wheel
pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts
pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq pyyaml

## Install R
apt-get install -y r-base cmake

# configure emulandice
$FACTSROOT/modules/emulandice/emulandice_config.sh

# set up radical pilot sandbox
# 
# If you do not create symlinks here, then running as localhost,
# the radical.pilot.sandbox directory (and the potentially very large files it contains) will reside in the
# user's home directory.

if [[ ! -z $RADICALPILOTSANDBOX ]]; then
    if [ ! -d $RADICALPILOTSANDBOX ]; then
        mkdir $RADICALPILOTSANDBOX
    fi

    if [ ! -d $HOME/radical.pilot.sandbox ]; then
        ln -s  $RADICALPILOTSANDBOX $HOME/radical.pilot.sandbox
    fi
fi

# do a test run

cd /opt/facts
python runFACTS.py experiments/dummy
