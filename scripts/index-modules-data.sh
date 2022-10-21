#!/bin/bash

for i in *.tgz
do
   echo
   echo $i
   echo -------------------------------
   tar tzf $i 2>&1 | grep -v 'Ignoring'
done
