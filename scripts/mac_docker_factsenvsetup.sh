#!/bin/bash

# This script will create a Docker environment that can run FACTS, starting with the command
#
# docker run --hostname=localhost --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch --runtime=runc -it ubuntu:focal

apt-get update -y
apt-get install -y sudo systemctl gnupg wget curl apt-transport-https
wget -qO - https://www.mongodb.org/static/pgp/server-6.0.asc | sudo apt-key add -
echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/6.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-6.0.list
apt-get update -y
apt-get install -y mongodb-org

## Install python
apt-get install -y python3-pip python3.8-venv git libnetcdf-dev python3-netcdf4
apt-get install -y ssh iputils-ping slurm

python3 -m venv ve3
. ve3/bin/activate
pip install --upgrade setuptools pip wheel
pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts
pip install numpy scipy netCDF4 pyyaml matplotlib h5py yq pyyaml

systemctl start mongod

sleep 5

# set up radical pilot sandbox

if [ ! -d /scratch/radical.pilot.sandbox ]; then
    mkdir  /scratch/radical.pilot.sandbox
fi

if [ ! -d $HOME/radical.pilot.sandbox ]; then
    ln -s  /scratch/radical.pilot.sandbox $HOME/radical.pilot.sandbox
fi

# do a test run

cd /opt/facts
python runFACTS.py experiments/dummy
