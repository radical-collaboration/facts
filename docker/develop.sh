#!/bin/bash

TAG="facts"

echo "Build docker container $TAG"
docker build --no-cache --target facts-core -t "$TAG" .

TAG="facts-jupyter"
echo "Build docker container $TAG"
docker build --no-cache --target facts-jupyter -t "$TAG" .
