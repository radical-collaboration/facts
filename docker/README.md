# Executing FACTS on a Docker container

To install FACTS through docker please follow the steps below:

1. from the root of FACTS repository:

   ```shell
   cd docker
   ```

2. Build the docker container image, naming it `facts`:

   ```shell
   sh develop.sh
   ```

3. Create a directory for the RADICAL Pilot sandbox:

   ```shell
   mkdir -p ~/tmp/radical.pilot.sandbox
   ```

4. Start a container from the `facts` image. Note: you will have to specify the directory where you cloned the FACTS repository. The following assume that the FACTS repository was cloned in `$HOME/facts`:

   ```shell
   docker run --hostname=localhost \
              --runtime=runc \
              --volume=$HOME/facts:/opt/facts \
              --volume=$HOME/tmp/radical.pilot.sandbox:/root/radical.pilot.sandbox \
              -it \
              facts
   ```

5. Confirm that FACTS work within the container:

   ```shell
   cd /opt/facts ; python3 runFACTS.py experiments/dummy
   ```