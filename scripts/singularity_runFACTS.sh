#!/bin/bash

FACTSDIR=$HOME/facts
LOGDIR=$FACTSDIR/log

if [ ! -d $LOGDIR ]; then
    mkdir $LOGDIR
fi

singularity instance start -w $HOME/singularity/rabbitmq rabbitmq
singularity instance start -w $HOME/singularity/mongo mongo

singularity run -w instance://rabbitmq > $LOGDIR/rabbitmq.`hostname`.out 2> $LOGDIR/rabbitmq.`hostname`.err &
singularity run -w instance://mongo > $LOGDIR/mongo.`hostname`.out 2> $LOGDIR/mongo.`hostname`.err &

sleep 30

j=0
for i in $@
do
    python $FACTSDIR/runFACTS.py $i > $LOGDIR/experiment.$j.`hostname`.out 2> $LOGDIR/experiment.$j.`hostname`.err
    j=$(($j + 1))
done
