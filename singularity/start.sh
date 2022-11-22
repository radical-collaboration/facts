#!/bin/sh

tmp="$(pwd)/rp_tmp/"
scratch="$(pwd)/scratch/"
container='rp_services'

mkdir -p $tmp
chmod 1777 $tmp

mongodb_conf=mongodb.conf
rabbitmq_conf=rabbitmq.conf

cp -f $mongodb_conf  $tmp/
cp -f $rabbitmq_conf $tmp/

# get unique port numbers for each user
uid=$(id -u)
mongodb_port=$((  10000 + (uid - 1000) * 3 + 0))
rabbitmq_nport=$((10000 + (uid - 1000) * 3 + 1))
rabbitmq_dport=$((10000 + (uid - 1000) * 3 + 2))
rabbitmq_nname="RMQ_$rabbitmq_nport"

# configure port numbers in mongodb rabbitmq config files
sed -i -e "s/###MONGODB_PORT###/$mongodb_port/g"      $tmp/$mongodb_conf
sed -i -e "s/###RABBITMQ_CONF###/$rabbitmq_conf/g"    $tmp/$rabbitmq_conf
sed -i -e "s/###RABBITMQ_NPORT###/$rabbitmq_nport/g"  $tmp/$rabbitmq_conf
sed -i -e "s/###RABBITMQ_DPORT###/$rabbitmq_dport/g"  $tmp/$rabbitmq_conf
sed -i -e "s/###RABBITMQ_NNAME###/$rabbitmq_nname/g"  $tmp/$rabbitmq_conf

# write user environment
cat <<EOT > $tmp/$container.env

export RADICAL_PILOT_DBURL=mongodb://localhost:$mongodb_port/rct

export RABBITMQ_CONFIG_FILE=/tmp/rabbitmq.conf
export RMQ_HOST=localhost
export RMQ_PORT=$rabbitmq_port
export RMQ_USER=''
export RMQ_PASS=''

EOT

# add RMQ settings to container env
cat $tmp/$rabbitmq_conf | sed -e 's/\(^.*=\)/export \1/g' >> $tmp/$container.env


# only build new container if it doesn't exist yet
if ! test -d "$scratch/$container"
then
    singularity build \
        --fakeroot \
        --force \
        --bind $tmp:/tmp/ \
        --sandbox $scratch/$container \
        $container.spec
fi

# start the container and run services
singularity instance start \
    --bind $tmp:/tmp \
    $scratch/$container rps

singularity run \
    --bind $tmp:/tmp \
    instance://rps

