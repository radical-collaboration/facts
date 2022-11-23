#!/bin/sh

scratch="$(pwd)/scratch/"
tmp="$(pwd)/rct_tmp/"

mkdir -p $scratch
mkdir -p $tmp
chmod 1777 $tmp

mongodb_conf=mongodb.conf

cp -f $mongodb_conf $tmp/

# get unique port numbers for each user
uid=$(id -u)
mongodb_port=$((  10000 + (uid - 1000)))

# configure port numbers in mongodb rabbitmq config files
sed -i -e "s/###MONGODB_PORT###/$mongodb_port/g"      $tmp/$mongodb_conf

# only build new container if it doesn't exist yet
if ! test -d "$scratch/rct_mongodb"
then
    singularity --bind $tmp:/tmp/ build \
        --fakeroot \
	--force \
        --bind $tmp:/tmp/ \
        --sandbox $scratch/rct_mongodb \
        rct_mongodb.spec
fi

if ! test -d "$scratch/rct_rabbitmq"
then
    singularity build \
        --sandbox $scratch/rct_rabbitmq docker://rabbitmq
fi

# make sure no old containers are running
singularity instance list | grep rct_mongodb  \
    && singularity instance stop rct_mongodb
singularity instance list | grep rct_rabbitmq \
    && singularity instance stop rct_rabbitmq

# start the container and run services
singularity instance start \
    --bind $tmp:/tmp \
    $scratch/rct_mongodb rct_mongodb

singularity instance start \
    --bind $tmp:/tmp \
    $scratch/rct_rabbitmq rct_rabbitmq


singularity run \
    --bind $tmp:/tmp \
    $scratch/rct_mongodb \
    > $tmp/rct_mongodb.log 2>&1 &

singularity run \
    -w $scratch/rct_rabbitmq \
    > $tmp/rct_rabbitmq.log 2>&1 &
echo "$!" > $tmp/rct_rabbitmq.pid

# grab rabbitmq connection details from logfile
echo "check for service endpoints"
start=$(date +%s)
while true
do
    now=$(date +%s)
    diff=$((now - start))
    test "$diff" -gt 10 && break
    port=$(cat $tmp/rct_rabbitmq.log \
             | grep 'started TCP listener on' \
             | awk -F : '{print $NF}' \
             | sed -e 's/.*//g')   #  remove trailing ansi escapes
    test -z "$port" || break
    sleep 1
done

echo "write rct_services.env"
cat <<EOT > rct_services.env

export RADICAL_PILOT_DBURL=mongodb://localhost:$mongodb_port/rct

export RMQ_HOST=localhost
export RMQ_PORT=$port
export RMQ_USER=
export RMQ_PASS=

EOT

