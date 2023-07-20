FROM ubuntu:focal
WORKDIR /opt/facts
ARG MONGODB_URL=mongodb://guest:guest@mongodb:27017/default

# Allows virtual enviorment to be activated
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# updates packages and installs python3 as well as git
RUN apt-get update
RUN apt-get install -y python3.9
RUN apt-get update && apt-get install -y git

# Install pip from standard ubuntu packages
RUN apt-get install -y python3-pip
RUN apt-get install  -y python3.8-venv

# Installs and starts mongodb
# NOTE: COMMENTED MONGODB INSTALL AT BOTTOM REQUIRES HANDLING OF USER INPUT
RUN apt-get install -y sudo gnupg curl apt-transport-https
RUN curl -fsSL https://pgp.mongodb.com/server-6.0.asc | \
   sudo gpg -o /usr/share/keyrings/mongodb-server-6.0.gpg \
   --dearmor
RUN echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/6.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-6.0.list
RUN apt-get update -y
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends mongodb-org

# Starts MongoDB
RUN /usr/bin/mongod --config /etc/mongod.conf --fork

# Creates and activates python3 virtual environment 
RUN python3 -m venv factsVe
RUN source factsVe/bin/activate

#Pip installs required FACTS packages
RUN pip install --upgrade setuptools pip wheel && \
pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts && \
pip install numpy==1.24.3 && \
pip install scipy==1.10.1 && \ 
pip install netCDF4==1.6.3 && \
pip install pyyaml==6.0 && \
pip install matplotlib==3.7.1 && \ 
pip install h5py==3.8.0 && \
pip install yq==3.2.2

# Installs R for ubuntu20.04 Focal
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends r-base cmake