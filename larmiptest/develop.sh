#!/bin/bash

TAG="larmiptest"

echo "Build docker container $TAG"
docker build -t "$TAG" .

# docker tag "$TAG":latest factscapstone/facts:"$TAG"
# docker push factscapstone/facts:"$TAG"