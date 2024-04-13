#!/bin/sh

find totaltest/to_total -type f -delete
cp totaltest/configTest.ocean.tlm.sterodynamics_globalsl.nc totaltest/to_total/global/
cp totaltest/configTest.ocean.tlm.sterodynamics_localsl.nc totaltest/to_total/local/

sed "s+homedir+$(echo $HOME)+g" global_kubernetes.yaml > temp_file.yaml && mv temp_file.yaml kubernetes.yaml

find *test -name "develop.sh" -exec sh {} \;

kubectl delete pod test
kubectl apply -f kubernetes.yaml
kubectl get pods
