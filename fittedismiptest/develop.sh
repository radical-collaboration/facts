#!/bin/bash

TAG="fittedismiptest"

echo "Build docker image $TAG"
docker build -t "$TAG" "$HOME"/Capstone/facts/"$TAG"

# echo "Start the container"
# docker run -it fittedismiptest

# echo "Push the image to docker registry"
