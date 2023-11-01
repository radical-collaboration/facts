#!/bin/bash
trun="${WALLTIME:-600}"
srun -p kopp_1 -c 1 --time=$trun --propagate=NONE python3 runFACTS.py $1