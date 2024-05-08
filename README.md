# Framework for Assessing Changes To Sea-level (FACTS)

## Capstone Branch

This branch contains work that was submitted in partial fulfillment of the requirements for senior design project.

### Overview

A minimum viable product of the alternative solution for the infrastructure of FACTS is presented in this branch. In its current state, FACTS utilizes RADICAL EnTK to run its experiments. The new proposed system utilizes Docker and Kubernetes to manage the runtime as a replacement to EnTK.

Individual modules are housed in separate stand-alone directories with a Dockerfile and execute.sh script. This ensures that each module can be built into a Docker image that can be run independently as a Docker container. A developer that wants to dockerize their module would need to specify the dependencies required in the Dockerfile and the files to be executed in the execute.sh script.

To manage the modules as specified in an experiment, Kubernetes is implemented in this solution. Kubernetes allows the parallelization of modules in addition to the management of resource allocation for the running containers. For a researcher to outline an experiment, they would need to provide a kubernetes.yaml file. This file contains information for any shared volumes that will store data as well as the scheduling of modules to be run. Kubernetes will then build the containers based on the images provided and start the runtime for the experiment. An advantage of Kubernetes is that it is compatible with local machines, high-performance computing machines, and the cloud.

### How to Run

This minimum viable product encapsulated the following modules in an experiment: fair_temperature, fittedISMIP, stereodynamics, & larmipAIS.

First, we must install Docker Desktop and enable Kubernetes from the settings within Docker Desktop. Once the Docker engine as well as Kubernetes is running, the images for the required modules need to be created. For example,

- Navigate to the ```fairtest``` directory and execute ```sh ./develop.sh``` to generate the image *fairtest*.

The images must be generated for each of the modules that are to be run.

*Note: if the ```find``` command in the ```global.sh``` file is uncommented, the images will be generated automatically.*

To specify the modules that are to be run, the ```modules.yaml``` file can be edited. The containers are specified toward the end of the file. The structure to add containers remains the same with the only changes needing to be made is to the name of the image the container uses. **If a full run is desired of all four modules, this file needs not be updated.**

In the ```facts``` directory, the user can now run ```sh ./global.sh``` in order to run the jobs for each module.

Upon execution, the script will generate containers needed and the modules specified will be run. By clicking the containers, the logs from the run can be seen. Using either the ```exec``` or ```files``` tab in Docker Desktop, the outputs and files generated can be seen.

### Action Items

- Dockerize remaining modules.
- To create more robust operation, the ```execute.sh``` files in each module directory need to be updated to use a specific number of locations and samples.

[![DOI](https://zenodo.org/badge/151614681.svg)](https://zenodo.org/badge/latestdoi/151614681)

The Framework for Assessing Changes To Sea-level (FACTS) is an open-source modular, scalable, and extensive framework for global mean, regional, and extreme sea level projection that is designed to support the characterization of ambiguity in sea-level projections. It is designed so users can easily explore deep uncertainty by investigating the implications on GMSL, RSL, and ESL of different choices for different processes. Its modularity allows components to be represented by either simple or complex model. Because it is built upon the Radical-PILOT computing stack, different modules can be dispatched for execution on resources appropriate to their computational complexity.

FACTS is being developed by the [Earth System Science & Policy Lab](https://www.earthscipol.net) and the [RADICAL Research Group](https://radical.rutgers.edu) at Rutgers University. FACTS is released under the MIT License.

See [fact-sealevel.readthedocs.io](https://fact-sealevel.readthedocs.io) for documentation.
