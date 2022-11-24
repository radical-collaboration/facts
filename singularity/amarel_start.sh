#!/bin/sh

module load singularity

scratch="/scratch/$(id -un)/singularity"
tmp="$scratch/tmp"

mkdir -p   $tmp/db
chmod 1777 $tmp

# only build new container if they don't exist yet
if ! test -d "$scratch/rct_mongodb"
then
    singularity build --sandbox $scratch/rct_mongodb docker://mongo
fi

if ! test -d "$scratch/rct_rabbitmq"
then
    singularity build --sandbox $scratch/rct_rabbitmq docker://rabbitmq
fi

# make sure no old containers are running
singularity instance list | grep rct_mongodb \
    && singularity instance stop rct_mongodb

singularity instance list | grep rct_rabbitmq \
    && singularity instance stop rct_rabbitmq

# start the container
singularity instance start \
    $scratch/rct_mongodb rct_mongodb

singularity instance start \
    $scratch/rct_rabbitmq rct_rabbitmq

# run services
singularity run \
    --bind $tmp:/tmp \
    --bind $tmp:/data \
    $scratch/rct_mongodb mongod \
    > $tmp/rct_mongodb.log 2>&1 &

singularity run \
    -w $scratch/rct_rabbitmq \
    > $tmp/rct_rabbitmq.log 2>&1 &

