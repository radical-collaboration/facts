#!/bin/bash

TAG="larmiptest"

echo "Build docker container $TAG"
docker build -t "$TAG" .

docker tag "$TAG":latest $(whoami)/"$TAG":latest
docker push $(whoami)/"$TAG":latest
