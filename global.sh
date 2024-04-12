#!/bin/sh

sed "s/whoami/$(whoami)/g; s+homedir+$(echo $HOME)+g" kubernetes.yaml > temp_file.yaml && mv temp_file.yaml kubernetes.yaml

find *test -name "develop.sh" -exec sh {} \;

kubectl delete pod test
kubectl apply -f kubernetes.yaml
