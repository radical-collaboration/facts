#!/bin/bash

TAG="totaltest"

echo "Build docker container $TAG"
docker build -t "$TAG" .