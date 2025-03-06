#!/bin/bash

TAG="facts"

echo "Build docker container $TAG"
docker build --target facts-core -t "$TAG" .

TAG="facts-jupyter"
echo "Build docker container $TAG"
docker build --target facts-jupyter -t "$TAG" .
