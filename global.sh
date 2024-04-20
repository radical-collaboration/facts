#!/bin/bash

# sed "s+homedir+$(echo $HOME)+g" global_kubernetes.yaml > temp_file.yaml && mv temp_file.yaml kubernetes.yaml

# find *test -name "develop.sh" -exec sh {} \;

echo "Initializing persistent volumes and claims..."
kubectl apply -f init.yaml
kubectl wait --timeout=-1s --for=condition=complete job.batch/init-job
echo

echo "Starting modules pod..."
kubectl apply -f modules_job.yaml
echo

# echo "Waiting for first pod to complete..."
# kubectl wait --timeout=-1s --for=condition=complete job.batch/modules-job
# echo

# echo "Starting totaling pod..."
# kubectl apply -f total_job.yaml
