#!/bin/sh

tmp=`pwd`/rp_tmp/
scratch=`pwd`/scratch/
cont='rp_services'

singularity instance list | grep rps \
    && singularity instance stop rps
echo

rm -rvf $tmp/mongodb* $tmp/rabbitmq* rp_services.env
echo

