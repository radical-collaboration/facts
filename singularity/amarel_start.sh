#!/bin/sh

module load singularity

node=$(hostname | cut -f 1 -d .)
scratch="/scratch/$(id -un)/singularity"
tmp="$scratch/tmp_$node"

mkdir -p   $tmp/db
chmod 1777 $tmp

# only build new container if they don't exist yet
test -d "$scratch/rct_mongodb" \
  || singularity build --sandbox $scratch/rct_mongodb docker://mongo

test -d "$scratch/rct_rabbitmq" \
 || singularity build --sandbox $scratch/rct_rabbitmq docker://rabbitmq

# run services
singularity run \
    --bind $tmp:/tmp \
    --bind $tmp:/data \
    $scratch/rct_mongodb mongod \
    > $tmp/rct_mongodb.log 2>&1 &

singularity run \
    -w $scratch/rct_rabbitmq \
    > $tmp/rct_rabbitmq.log 2>&1 &

