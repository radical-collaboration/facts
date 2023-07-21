To install FACTS through docker please follow the steps below:

1. cd into the "docker" directory
2. Build the docker container image, naming it `facts`:
   >> sh develop.sh
3. Start a container from the `facts` image, assuming that the FACTS repository was cloned in `$HOME/facts`:
   >> docker run --hostname=localhost --runtime=runc -it  --volume=$HOME/facts:/opt/facts facts
4. Confirm that FACTS work within the container:
   >> cd /opt/facts ; python3 runFACTS.py experiments/dummy
