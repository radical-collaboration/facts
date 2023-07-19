#!/bin/bash

BASE_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )/../" &> /dev/null && pwd 2> /dev/null; )"

TAG="factsealevel/facts$TAG"

if [[ -z $PLATFORM ]]; then
    PLATFORM="linux/amd64"
fi

echo "Build docker container $TAG (name: ${NAME:-n/a})"
docker build \
    -t "$TAG" \
    --build-arg TUTORIAL_NAME="$NAME" \
    --build-arg BUILDPLATFORM="$PLATFORM" \
    -f "Dockerfile" \
    "src/"