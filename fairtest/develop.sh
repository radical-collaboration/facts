#!/bin/bash

TAG="fairtest"

echo "Build docker container $TAG"
docker build -t "$TAG" "$HOME"/Capstone/facts/"$TAG"

docker tag "$TAG":latest $(whoami)/"$TAG":latest
docker push $(whoami)/"$TAG":latest
