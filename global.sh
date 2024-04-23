#!/bin/bash

# sed "s+homedir+$(echo $HOME)+g" global_kubernetes.yaml > temp_file.yaml && mv temp_file.yaml kubernetes.yaml

find *test -name "develop.sh" -exec sh {} \;

echo "Initializing persistent volumes and claims..."
kubectl apply -f init.yaml
kubectl wait --timeout=-1s --for=condition=complete job.batch/init-job
echo

echo "Running fair module pod..."
kubectl apply -f fair_job.yaml
kubectl wait --timeout=-1s --for=condition=complete job.batch/fair-job
echo

echo "Running modules pod..."
kubectl apply -f modules_job.yaml
kubectl wait --timeout=-1s --for=condition=complete job.batch/modules-job
echo

echo "Running totaling pod..."
kubectl apply -f total_job.yaml
