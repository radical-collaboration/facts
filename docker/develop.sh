#!/bin/bash

TAG="facts"

echo "Build docker container $TAG"
docker build -t "$TAG" .