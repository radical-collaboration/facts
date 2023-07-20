To install FACTS through docker please follow the steps below:
NOTE: STAY IN /facts/docker while doing the install
1. cd into the "docker" directory
2. >> sh develop.sh (Builds the docker container as "facts")
3. To drop into the container (assuming facts is installed in $HOME/facts and $HOME/tmp exists)
   >> docker run --hostname=localhost --runtime=runc -it  --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch facts 
4. start mongodb from within the container:
   >> /usr/bin/mongod --config /etc/mongod.conf --fork
5. >> python3 runFACTS.py experiments/dummy (to confirm RCT/Mongodb install)
