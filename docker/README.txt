To install FACTS through docker please follow the steps below:

1. cd into the "docker" directory
2. >> sh develop.sh (Builds the docker container as "facts")
3. To drop into the container (assuming facts is installed in $HOME/facts and $HOME/tmp exists)
   >> docker run --hostname=localhost --runtime=runc -it  --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp:/scratch facts 
4. >> cd /opt/facts ; python3 runFACTS.py experiments/dummy (to confirm RCT/Mongodb install)
