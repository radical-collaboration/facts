# Framework for Assessing Changes To Sea-level (FACTS)

## Capstone Branch

This branch contains work that was submitted in partial fulfillment of the requirements for senior design project. 

### Overview

A minimum viable product of the alternative solution for the infrastructure of FACTS is presented in this branch. In its current state, FACTS utilizes RADICAL EnTK to run its experiments. The new proposed system utilizes Docker and Kubernetes to manage the runtime as a replacement to EnTK. 

Individual modules are housed in separate stand-alone directories with a Dockerfile and execute.sh script. This ensures that each module can be built into a Docker image that can be run independently as a Docker container. A developer that wants to dockerize their module would need to specify the dependencies required in the Dockerfile and the files to be executed in the execute.sh script. 

To manage the modules as specified in an experiment, Kubernetes is implemented in this solution. Kubernetes allows the parallelization of modules in addition to the management of resource allocation for the running containers. For a researcher to outline an experiment, they would need to provide a kubernetes.yaml file. This file contains information for any shared volumes that will store data as well as the scheduling of modules to be run. Kubernetes will then build the containers based on the images provided and start the runtime for the experiment. An advantage of Kubernetes is that it is compatible with local machines, high-performance computing machines, and the cloud. 

### How to Run
This minimum viable product encapsulated the following modules in an experiment: fair temperature, fitted ISMIP, stereodynamics, & larmip icesheet. 

- Install Docker desktop and Kubernetes. 
- Specify and update the path to the FACTS module data in the kubernetes.yaml file.
- Execute ```kubectl delete pod test``` to ensure that previous runs do not interfere.
- Execute ```kubectl apply -f kubernetes.yaml`` to start the experiment.
- Execute ```kubectl get pods``` to view the current running pods

### Action Items
- Enable the ability to specify global variables in the kubneretes.yaml file. 
- Dockerize remaining modules. 




[![DOI](https://zenodo.org/badge/151614681.svg)](https://zenodo.org/badge/latestdoi/151614681)

The Framework for Assessing Changes To Sea-level (FACTS) is an open-source modular, scalable, and extensive framework for global mean, regional, and extreme sea level projection that is designed to support the characterization of ambiguity in sea-level projections. It is designed so users can easily explore deep uncertainty by investigating the implications on GMSL, RSL, and ESL of different choices for different processes. Its modularity allows components to be represented by either simple or complex model. Because it is built upon the Radical-PILOT computing stack, different modules can be dispatched for execution on resources appropriate to their computational complexity.

FACTS is being developed by the [Earth System Science & Policy Lab](https://www.earthscipol.net) and the [RADICAL Research Group](https://radical.rutgers.edu) at Rutgers University. FACTS is released under the MIT License.

See [fact-sealevel.readthedocs.io](https://fact-sealevel.readthedocs.io) for documentation.
