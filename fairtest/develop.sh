#!/bin/bash

TAG="fairtest"

echo "Build docker container $TAG"
docker build -t "$TAG" .
