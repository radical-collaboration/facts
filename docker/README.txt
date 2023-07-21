To install FACTS through docker please follow the steps below:

1. cd into the "docker" directory
2. Build the docker container image, naming it `facts`:
   >> sh develop.sh
3. Create a directory for the RADICAL Pilot sandbox:
   >> mkdir -p ~/tmp/radical.pilot.sandbox
4. Start a container from the `facts` image, assuming that the FACTS repository was cloned in `$HOME/facts`:
   >> docker run --hostname=localhost --runtime=runc -it  --volume=$HOME/facts:/opt/facts --volume=$HOME/tmp/radical.pilot.sandbox:/root/radical.pilot.sandbox -w /opt/facts facts
5. Confirm that FACTS work within the container:
   >> python3 runFACTS.py experiments/dummy
