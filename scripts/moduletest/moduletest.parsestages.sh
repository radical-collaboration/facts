#!/bin/bash

STAGE_SCRIPT=( )
j=0
for i in ${STAGES[@]}
do
    j=$(( j+1 ))
    STAGE_SCRIPT[j]=`yq --raw-output  .\"$i\".task1.script ../pipeline.yml`
done
