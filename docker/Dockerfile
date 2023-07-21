FROM ubuntu:focal

ENV RADICAL_PILOT_DBURL="mongodb://127.0.0.1:27017/facts"
ENV RADICALPILOTSANDBOX="/scratch/radical.pilot.sandbox"
ENV DEBIAN_FRONTEND="noninteractive"

# allows virtualenv
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
        apt-transport-https

# Installs and starts mongodb
RUN curl -fsSL https://pgp.mongodb.com/server-6.0.asc | \
   sudo gpg -o /usr/share/keyrings/mongodb-server-6.0.gpg \
   --dearmor
RUN echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/6.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-6.0.list
RUN apt-get update -y
RUN apt-get install -y --no-install-recommends mongodb-org

# Creates and activates python3 virtual environment
RUN python3 -m venv factsVe &&\
    source factsVe/bin/activate

# Pip installs required FACTS packages
RUN pip install --upgrade \
        setuptools pip wheel \
        git+https://github.com/radical-cybertools/radical.entk@projects/facts \
        pyyaml==6.0

# Installs R for ubuntu20.04 Focal
RUN apt-get update
RUN apt-get install -y --no-install-recommends r-base cmake

# make directory for radical pilot sandbox
RUN mkdir -p ~/radical.pilot.sandbox

# Return a shell to the user, starting MDB
RUN echo "/usr/bin/mongod --config /etc/mongod.conf --fork" > ~/.bashrc
CMD /bin/bash
