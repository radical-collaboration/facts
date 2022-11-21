#!/bin/sh

tmp=`pwd`/my_tmp/
scratch=`pwd`/scratch/
cont='rp_services'

singularity instance stop rps
rm -f $tmp/mongodb.pid

