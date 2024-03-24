#!/bin/bash

TAG="larmiptest"

echo "Build docker container $TAG"
docker build -t "$TAG" .