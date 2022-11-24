#!/bin/bash

#SBATCH --partition=kopp_1
#SBATCH --nodes=3
#SBATCH --cpus-per-task=1  
#SBATCH --time=00:30:00 
#SBATCH --output=slurm.out   
#SBATCH --error=slurm.out

queue='kopp_1'
nnodes=3

cpn=$(cat /proc/cpuinfo | grep -e '^processor\>' | wc -l)
cores=$((nnodes * cpn))
hostname=$(hostname -f)

./amarel_start.sh
sleep 10
ps -ef | grep am3465
netstat -plunt


cd ~/facts/
. ve3/bin/activate

cat <<EOT > ./test/resource_amarel.yml

rabbitmq:
  hostname: "$hostname"

mongodb:
  hostname: "$hostname"

resource-desc:
  resource: 'rutgers.amarel'
  walltime: 30
  cpus: $cores
  gpus: 0
  queue: '$queue'
  project: ''
  schema: 'interactive'

EOT


python ./runFACTS.py --resourcedir test test

