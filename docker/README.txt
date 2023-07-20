To install FACTS through docker please follow the steps below:
NOTE: STAY IN /facts/docker while doing the install
1. cd into the "docker" directory
2. >> docker build . (runs the Dockerfile)
3. >> sh develop.sh (Builds the docker container as "facts")
4. To drop into the container 
5. >> docker run --hostname=localhost --runtime=runc -it -v $(pwd)/../opt/facts facts
6. start mongodb: /usr/bin/mongod --config /etc/mongod.conf --fork
7. >> python3 runFACTS.py experiments/dummy (to confirm RCT/Mongodb install)