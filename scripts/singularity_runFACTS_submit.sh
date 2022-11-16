#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=runFACTS
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --time=96:00:00
#SBATCH --output=out_runFACTS.txt
#SBATCH --error=error_runFACTS.txt
#SBATCH --export=ALL

cd $HOME/facts
source $HOME/.bashrc
conda activate facts
module load singularity

FACTSDIR=$HOME/facts
LOGDIR=$FACTSDIR/log

if [ ! -d $LOGDIR ]; then
    mkdir $LOGDIR
fi

singularity run -w $HOME/singularity/rabbitmq > $LOGDIR/rabbitmq.`hostname`.out 2> $LOGDIR/rabbitmq.`hostname`.err &
singularity run -w $HOME/singularity/mongo > $LOGDIR/mongo.`hostname`.out 2> $LOGDIR/mongo.`hostname`.err &

j=0
for i in $@
do
    python $FACTSDIR/runFACTS.py $i > $LOGDIR/experiment.$j.`hostname`.out 2> $LOGDIR/experiment.$j.`hostname`.err    
    j=$(($j + 1))
done
