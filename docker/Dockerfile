FROM ubuntu:focal

ENV RADICAL_PILOT_DBURL="mongodb://127.0.0.1:27017/facts"
ENV RADICALPILOTSANDBOX="/scratch/radical.pilot.sandbox"

# allows virtualenv
SHELL [ "/bin/bash", "-c" ]

# updates packages and installs python3 as well as git
RUN apt-get update
RUN apt-get install -y python3.9
RUN apt-get update && apt-get install -y git

# Install pip from standard ubuntu packages
RUN apt-get install -y python3-pip
RUN apt-get install  -y python3.8-venv
RUN apt-get install -y libhdf5-dev libnetcdf-dev

# Installs and starts mongodb
# NOTE: COMMENTED MONGODB INSTALL AT BOTTOM REQUIRES HANDLING OF USER INPUT
RUN apt-get install -y sudo gnupg curl apt-transport-https
RUN curl -fsSL https://pgp.mongodb.com/server-6.0.asc | \
   sudo gpg -o /usr/share/keyrings/mongodb-server-6.0.gpg \
   --dearmor
RUN echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb-server-6.0.gpg ] https://repo.mongodb.org/apt/ubuntu focal/mongodb-org/6.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-6.0.list
RUN apt-get update -y
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends mongodb-org


# Creates and activates python3 virtual environment 
RUN python3 -m venv factsVe
RUN source factsVe/bin/activate

#Pip installs required FACTS packages
RUN pip install --upgrade setuptools pip wheel && \
pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts && \
pip install pyyaml==6.0

# Installs R for ubuntu20.04 Focal
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends r-base cmake

RUN if [[ ! -z $RADICALPILOTSANDBOX ]]; then  \
    if [ ! -d $RADICALPILOTSANDBOX ]; then  \
        mkdir $RADICALPILOTSANDBOX ; \
    fi ; \
    if [ ! -d $HOME/radical.pilot.sandbox ]; then \
        ln -s  $RADICALPILOTSANDBOX $HOME/radical.pilot.sandbox ; \
    fi ; \
fi

RUN echo "/usr/bin/mongod --config /etc/mongod.conf --fork" > ~/.bashrc

CMD /bin/bash
