#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TAG="facts"

echo "Build docker container $TAG"
docker build --no-cache --target facts-core -t "$TAG" -f "${REPO_ROOT}/docker/Dockerfile" "${REPO_ROOT}"

TAG="facts-jupyter"
echo "Build docker container $TAG"
docker build --no-cache --target facts-jupyter -t "$TAG" -f "${REPO_ROOT}/docker/Dockerfile" "${REPO_ROOT}"
