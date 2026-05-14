#!/usr/bin/env bash
set -e

image_shared="/opt/facts-image/modules/emulandice/shared"
runtime_shared="/opt/facts/modules/emulandice/shared"
bundle="emulandice_bundled_dependencies.tgz"

if [[ -f "${image_shared}/${bundle}" && -d "${runtime_shared}" && ! -f "${runtime_shared}/${bundle}" ]]; then
    cp -a "${image_shared}/${bundle}" "${runtime_shared}/${bundle}"
fi

# Copy modules-data tarballs baked into the image into the (writable)
# runtime mount, but only the ones that aren't already present.
image_data="/opt/facts-image/modules-data"
runtime_data="/opt/facts/modules-data"
if [[ -d "${image_data}" && -d "${runtime_data}" ]]; then
    shopt -s nullglob
    for src in "${image_data}"/*.tgz; do
        name="$(basename "$src")"
        if [[ ! -e "${runtime_data}/${name}" ]]; then
            cp -a "$src" "${runtime_data}/${name}"
        fi
    done
    shopt -u nullglob
fi

exec "$@"
