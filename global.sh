#!/bin/sh

# sed "s+homedir+$(echo $HOME)+g" global_kubernetes.yaml > temp_file.yaml && mv temp_file.yaml kubernetes.yaml

# find *test -name "develop.sh" -exec sh {} \;

echo "Deleting files in persistent volume..."
kubectl apply -f delete_files_job.yaml
kubectl wait --timeout=-1s --for=condition=complete job.batch/delete-files-job
echo

echo "Starting modules pod..."
kubectl apply -f modules_job.yaml
echo

echo "Waiting for first pod to complete..."
kubectl wait --timeout=-1s --for=condition=complete job.batch/modules-job
echo

echo "Starting totaling pod..."
kubectl apply -f total_job.yaml
