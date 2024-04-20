FROM ubuntu:focal

# set env variables for RCT and apt
ENV DEBIAN_FRONTEND="noninteractive"

# set shell for RUN instruction
SHELL [ "/bin/bash", "-c" ]

# installs dependences
RUN apt-get update &&\
    apt-get install -y \
        python3.9 \
        python3-pip \
        python3.8-venv \
        git \
        sudo \
        libhdf5-dev \
        libnetcdf-dev \
        gnupg \
        curl \
        neovim \
        apt-transport-https

# Creates and activates python3 virtual environment
RUN python3 -m venv --system-site-packages factsVe &&\
    source factsVe/bin/activate

# installs required FACTS packages
RUN pip install --upgrade \
        setuptools pip wheel \
        radical.entk \
        pyyaml==6.0 \
        memory_profiler \
        numpy \
        scipy \
        netCDF4 \
        matplotlib \
        pandas \
        xarray \
        dask \
        h5py \
        fair==1.6.4

# installs R for ubuntu2 0.04 Focal
RUN apt-get update &&\
    apt-get install -y --no-install-recommends r-base \
        cmake \
        libopenblas-dev \
        gfortran

# make directory for radical pilot sandbox
RUN mkdir -p ~/radical.pilot.sandbox

COPY ./ /opt/facts
WORKDIR /opt/facts

# Activate venv
RUN echo -e ". /factsVe/bin/activate" > ~/.bashrc

# return a shell to the user
CMD /bin/bash
