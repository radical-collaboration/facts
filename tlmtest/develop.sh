#!/bin/bash

TAG="tlmtest"

echo "Build docker container $TAG"
docker build -t "$TAG" .