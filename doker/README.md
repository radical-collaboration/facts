# Docker container with FACTS

## 1. Get container image

(*) The FACTS container is based on [jupyter/minimal-notebook](https://github.com/jupyter/docker-stacks) image.

**To have the FACTS container image locally it should be either BUILT or PULLED from DockerHub.**

### 1.A. Build container

```shell
./docker/build.sh
```

### 1.B. Pull container

```shell
docker pull facts/facts:latest
```

## 2. Run container image

### 2.A. Run `docker-compose`

It starts the `facts` container with the auxiliary service MongoDB, which is used by the RADICAL-Cybertools components as part of a communication layer.

```shell
cd docker
docker compose up -d
docker compose logs -f facts
# stop containers
#    docker compose stop
# remove containers
#    docker compose rm -f
```

### 2.B. Run container image with MongoDB service manually

These steps do the same as `docker-compose`, but all necessary commands are executed manually.

Docker network to communicate with service(s):

```shell
docker network create rct-network
```

Launch MongoDB service:

```shell
docker run -d --hostname mongodb --name rct-mongodb -p 27017:27017 \
           -e MONGO_INITDB_ROOT_USERNAME=root_user \
           -e MONGO_INITDB_ROOT_PASSWORD=root_pass \
           -e MONGO_INITDB_USERNAME=guest \
           -e MONGO_INITDB_PASSWORD=guest \
           -e MONGO_INITDB_DATABASE=default \
           --network rct-network mongo:4.4
```
```shell
docker exec rct-mongodb bash -c \
  "mongo --authenticationDatabase admin -u root_user -p root_pass default \
   --eval \"db.createUser({user: 'guest', pwd: 'guest', \
                           roles: [{role: 'readWrite', db: 'default'}]});\""
```

Run container with network:

```shell
docker run --rm -it -p 8888:8888 --network rct-network facts/facts
```

Stop services after work is done:

```shell
# stop container(s)
docker stop rct-mongodb
# stop and remove container(s)
#    docker rm -f rct-mongodb
```
