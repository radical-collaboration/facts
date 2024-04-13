#!/bin/bash

TAG="totaltest"

echo "Build docker container $TAG"
docker build -t "$TAG" "$HOME"/Capstone/facts/"$TAG"

docker tag "$TAG":latest $(whoami)/facts:"$TAG"
docker push $(whoami)/facts:"$TAG"
