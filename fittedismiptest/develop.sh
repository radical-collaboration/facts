#!/bin/bash

TAG="fittedismiptest"

echo "Build docker container $TAG"
docker build -t "$TAG" .
