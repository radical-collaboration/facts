#!/bin/bash

TAG="fairtest"

echo "Build docker container $TAG"
docker build -t "$TAG" "$HOME"/Capstone/facts/"$TAG"

docker tag "$TAG":latest factscapstone/facts:"$TAG"
docker push factscapstone/facts:"$TAG"
