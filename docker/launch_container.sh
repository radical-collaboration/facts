#!/usr/bin/env bash

# -Usage----------------------------------------------------------------
# This script creates and/or launches a FACTS Docker image/container.
#
# Run using:
#   bash launch_container.sh 
#
# Before running,
# review and update the user configuration in STEP 0,
# ==> especially:
#   MODE, IMAGE , container_name, CPU, memory,
#   facts_modules_data
#
# All paths are resolved relative to the FACTS repo root (the parent
# directory of this script).
# -------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"


#-STEP 0-------------------------------------------------------------
#         User configuration.
#--------------------------------------------------------------------

# Select ONE mode ONLY:
MODE="${MODE:-full}"      # build FACTS image & launch the container. (First time using) 
# MODE="${MODE:-run}"     # launch the container using an existing image

# Main Docker image name (e.g.: ssisls)
IMAGE="${IMAGE:-ssisls}"

# Container name & details
container_name="${IMAGE}_$(date +%Y%m%d_%H%M)"
CPU="${CPU:-8}"                                   # CPU (in terminal , linux: nproc    , mac:`sysctl hw.ncpu`)
memory="${memory:-12g}"                           # RAM (in terminal , linux: free -h  , mac:`system_profiler SPHardwareDataType | grep "Memory:"`)  

# Path to FACTS modules-data directory (relative paths resolve against REPO_ROOT)
facts_modules_data="${facts_modules_data:-modules-data}"
# facts_modules_data="${facts_modules_data:-/Users/uname/Desktop/FACTS_dev/modules-data}"    # Use for alternate location for data



#- End of user configuration-----------------------------------------
#  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
#- Below is ONLY for advanced users! 
# -------------------------------------------------------------------

# Sandbox options
sandbox_path="${sandbox_path:-_scratch/radical.pilot.sandbox}"
sandbox="${sandbox:-sandbox_path}"   
# sandbox="${sandbox:-tmp}"  # tmp | docker_volume_sandbox

# Bake modules-data tarballs into the image at build time:
#   none   = skip (default; fastest build)
#   global = download the global-only URL list  # Using this is slow 
#   all    = download the full URL list.        # Using this is slow 
MODULES_DATA="${MODULES_DATA:-none}"

#-Helpers---------------------
# 
# ----------------------------
die() { echo "ERROR: $*" >&2; exit 1; }

validate_config() {
  [[ "$memory" =~ ^[0-9]+[bkmgBKMG]?$ ]] || die "memory='$memory' is invalid. Use a Docker size like 5g, 30g."
}

verify_repo_root() {
  [[ -d "${REPO_ROOT}/modules" && -d "${REPO_ROOT}/docker" ]] \
    || die "REPO_ROOT='${REPO_ROOT}' does not look like a FACTS repo (missing modules/ or docker/)."
}

resolve_path() {
  # Echo $1 unchanged if absolute; otherwise prepend REPO_ROOT.
  if [[ "$1" = /* ]]; then printf '%s\n' "$1"; else printf '%s/%s\n' "$REPO_ROOT" "$1"; fi
}

banner() {
cat <<'EOF'
#############################################
#                                           #
#   Welcome to the FACTS docker container   #
#                                           #
#############################################
EOF
}


#-STEP 1-------------------------------------------------------------
#         Check Docker daemon status.
#--------------------------------------------------------------------
require_docker() {
  docker info >/dev/null 2>&1 || die "Docker daemon is not running... launch it."
  echo "Docker daemon is running."
}


#-STEP 2-------------------------------------------------------------
#         Build docker Image. 
#--------------------------------------------------------------------
build_images() {
  case "$MODULES_DATA" in
    none|global|all) ;;
    *) die "MODULES_DATA='$MODULES_DATA' invalid. Use one of: none, global, all." ;;
  esac

  echo "Building docker Image $IMAGE  (MODULES_DATA=$MODULES_DATA)"
  docker build --no-cache --target facts-core \
    --build-arg "MODULES_DATA=${MODULES_DATA}" \
    -t "$IMAGE" -f "${REPO_ROOT}/docker/Dockerfile" "${REPO_ROOT}"

  if [[ -n "${IMAGE1:-}" ]]; then
    echo "Building docker Image $IMAGE1  (MODULES_DATA=$MODULES_DATA)"
    docker build --no-cache --target facts-jupyter \
      --build-arg "MODULES_DATA=${MODULES_DATA}" \
      -t "$IMAGE1" -f "${REPO_ROOT}/docker/Dockerfile" "${REPO_ROOT}"
  fi
}


#-STEP 3-------------------------------------------------------------
#         Create a docker volume for radical sandbox.
#--------------------------------------------------------------------
ensure_sandbox_volume() {
  [[ "$sandbox" == "docker_volume_sandbox" ]] || return 0

  if docker volume inspect facts_sandbox >/dev/null 2>&1; then
    echo "Volume 'facts_sandbox' already exists, skipping creation."
  else
    docker volume create facts_sandbox >/dev/null || die "Failed to create volume 'facts_sandbox'."
    echo "Created volume 'facts_sandbox'."
  fi
}



#-STEP 4-------------------------------------------------------------
#         Launch a docker Container.
#--------------------------------------------------------------------
run_container() {
  local sandbox_mount
  local modules_data_abs
  modules_data_abs="$(resolve_path "$facts_modules_data")"

  case "$sandbox" in
    docker_volume_sandbox)
      sandbox_mount="--volume=facts_sandbox:/home/jovyan/radical.pilot.sandbox"
      ;;
    tmp)
      mkdir -p "${REPO_ROOT}/tmp/radical.pilot.sandbox"
      sandbox_mount="--volume=${REPO_ROOT}/tmp/radical.pilot.sandbox:/home/jovyan/radical.pilot.sandbox"
      ;;
    sandbox_path)
      sandbox_path="$(resolve_path "$sandbox_path")"
      mkdir -p "$sandbox_path"
      sandbox_mount="--volume=${sandbox_path}:/home/jovyan/radical.pilot.sandbox"
      ;;
    *) die "Unknown sandbox='$sandbox'";;
  esac

  # Common args (arrays avoid quoting bugs)
  local -a run_args=(
    -it --init
    --name "$container_name"
    --cpus "$CPU"
    --memory "$memory"
    --memory-swap "$memory"
    -e HDF5_USE_FILE_LOCKING=FALSE
    --volume "${REPO_ROOT}:/opt/facts"
    --volume "${modules_data_abs}:/opt/facts/modules-data"
    -w /opt/facts
  )

  docker run "${run_args[@]}" "$sandbox_mount" "$IMAGE" /bin/bash
}


#-STEP 5-------------------------------------------------------------
#         MAIN.
#--------------------------------------------------------------------
printf '\n\n'
require_docker
verify_repo_root
validate_config
printf '\n\n'

case "$MODE" in
  full)
    build_images
    printf '\n\n'
    ;;
  run) ;;
  *)
    die "Unknown MODE='$MODE' (expected 'full' or 'run')"
    ;;
esac

ensure_sandbox_volume
printf '\n\n'

banner
printf '\n\n'
run_container
