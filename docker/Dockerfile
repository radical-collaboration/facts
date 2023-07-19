
# https://docs.docker.com/engine/reference/builder/#automatic-platform-args-in-the-global-scope
ARG BUILDPLATFORM=linux/amd64
ARG PY_VER=3.9

FROM --platform=${BUILDPLATFORM} jupyter/minimal-notebook:python-${PY_VER}

#ARG TUTORIAL_NAME=default
ARG MONGODB_URL=mongodb://guest:guest@mongodb:27017/default

ENV JUPYTER_ENABLE_LAB=yes
ENV RADICAL_PILOT_DBURL=${MONGODB_URL}

USER root

# system packages
RUN apt-get update  -y && \
    apt-get install -y curl gnupg sudo

# add user to sudoers
RUN adduser ${NB_USER} sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' > /etc/sudoers.d/no_pass

# config kernel(s)
RUN mamba update -n base -y -c conda-forge conda \
 && mamba clean -a -f -y \
 && printf '#!/bin/bash \n\
\n\
eval "$(conda shell.posix hook)" \n\
exec "$@" \n\
' > "${CONDA_DIR}/bin/kernelstart.sh" \
 && chmod +x "${CONDA_DIR}/bin/kernelstart.sh" \
 && sed -i 's/\[/\[\"kernelstart.sh\",/g' \
    "${CONDA_DIR}/share/jupyter/kernels/python3/kernel.json" \
 && fix-permissions "${CONDA_DIR}" \
 && fix-permissions "/home/${NB_USER}"

# Installs R for ubuntu20.04 Focal
#RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends r-base cmake
RUN apt-get update  -y && \
    apt-get install -y r-base cmake

USER ${NB_UID}

# install RCT components (base conda env)
RUN mamba install -y -n base \
    radical.pilot \
    radical.entk \
    radical.analytics \
 && mamba clean -a -f -y

# Pip installs required FACTS packages
# NOTE: CANNOT GET REQUIRMENTS FILE TO WORK!!!
RUN pip install --upgrade setuptools pip wheel && \
pip install git+https://github.com/radical-cybertools/radical.entk@projects/facts && \
pip install numpy==1.24.3 && \
pip install scipy==1.10.1  && \
pip install netCDF4==1.6.3 && \
pip install pyyaml==6.0 && \
pip install matplotlib==3.7.1 && \
pip install h5py==3.8.0 && \
pip install yq==3.2.2 

#Pull the facts repository from github
RUN git clone https://github.com/radical-collaboration/facts.git

# tutorials
#WORKDIR /facts/modules/emulandice/shared/

#RUN source ./facts/modules/emulandice/shared/emulandice_environment.sh && \ 
#RUN source ./facts/modules/emulandice/shared/emulandice_build.sh && \
RUN cd ./facts/modules/emulandice/shared/ && sudo ./emulandice_bundle_dependencies.sh
#COPY --chown=${NB_UID}:${NB_GID} ./${TUTORIAL_NAME} ./

RUN if [[ -f ./setup.sh ]] ; then \
       bash setup.sh && rm setup.sh ; \
    fi

RUN if [[ -f ./environment.yml ]] ; then \
       mamba env update -n base -f ./environment.yml && rm environment.yml ; \
    fi