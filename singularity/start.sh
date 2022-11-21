#!/bin/sh

tmp="$(pwd)/rp_tmp/"
scratch="$(pwd)/scratch/"
container='rp_services'

mkdir -p $tmp
chmod 1777 $tmp

mongodb_conf=mongodb.conf
rabbitmq_conf=rabbitmq.conf

cp -f $mongodb_conf  $tmp
cp -f $rabbitmq_conf $tmp

# get unique port numbers for each user
uid=$(id -u)
mongodb_port=$(( 10000 + uid * 2 + 0))
rabbitmq_port=$((10000 + uid * 2 + 1))

# configure port numbers in mongodb rabbitmq config files
sed -i -e "s/###MONGODB_PORT###/$mongodb_port/g"   $tmp/$mongodb_conf
sed -i -e "s/###RABBITMQ_PORT###/$rabbitmq_port/g" $tmp/$rabbitmq_conf
sed -i -e "s/###RABBITMQ_CONF###/$rabbitmq_port/g" $tmp/$rabbitmq_conf

# write user environment
cat <<EOT > $tmp/$container.env

export RADICAL_PILOT_DBURL=mongodb://localhost:$mongodb_port/rct

export RABBITMQ_CONFIG_FILE=/tmp/rabbitmq_conf
export RMQ_HOST=localhost
export RMQ_PORT=$rabbitmq_port
export RMQ_USER=''
export RMQ_PASS=''

EOT


# only build new container if it doesn't exist yet
if ! test -d "$scratch/$container"
then
    singularity build \
        --fakeroot \
        --force \
        --bind "`pwd`/my_tmp/:/tmp/" \
        --sandbox $scratch/$container \
        $container.spec
fi

singularity instance start \
    --bind `pwd`/my_tmp:/tmp \
    $scratch/$container rps

singularity run \
    --bind `pwd`/my_tmp:/tmp \
    instance://rps

