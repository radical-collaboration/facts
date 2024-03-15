#!/bin/sh

# ------------------------------------------------------------------------------

export RP_TASK_ID="task.000000"
export RP_TASK_NAME="task.000000,task1,stage.0000,preprocess,pipeline.0000,test.GrIS1f.FittedISMIP.GrIS"
export RP_PILOT_ID="pilot.0000"
export RP_SESSION_ID="re.session.c27d8e62-e314-11ee-a7e7-0242ac110002"
export RP_RESOURCE="local.localhost"
export RP_RESOURCE_SANDBOX="/root/radical.pilot.sandbox"
export RP_SESSION_SANDBOX="$RP_RESOURCE_SANDBOX/$RP_SESSION_ID/"
export RP_PILOT_SANDBOX="$RP_SESSION_SANDBOX/pilot.0000/"
export RP_TASK_SANDBOX="$RP_PILOT_SANDBOX/task.000000"
export RP_REGISTRY_ADDRESS="tcp://172.17.0.2:10002"
export RP_CORES_PER_RANK=1
export RP_GPUS_PER_RANK=0
export RP_GTOD="$RP_PILOT_SANDBOX/gtod"
export RP_PROF="$RP_PILOT_SANDBOX/prof"
export RP_PROF_TGT="$RP_PILOT_SANDBOX/task.000000/task.000000.prof"

rp_error() {
    echo "$1 failed" 1>&2
    exit 1
}

# ------------------------------------------------------------------------------
# rank ID
export RP_RANKS=1
export RP_RANK=0

rp_sync_ranks() {
    sig=$1
    echo $RP_RANK >> $sig.sig
    while test $(cat $sig.sig | wc -l) -lt $RP_RANKS; do
        sleep 1
    done
}

# ------------------------------------------------------------------------------
$RP_PROF exec_start ""

# ------------------------------------------------------------------------------
# pre-exec commands
$RP_PROF exec_pre ""
. $RP_PILOT_SANDBOX/env/rp_named_env.rp.sh || true || rp_error pre_exec
tar -xvf ipccar6_climate_data.tgz 2> /dev/null; rm ipccar6_climate_data.tgz || rp_error pre_exec
pip install --upgrade pip; pip install numpy scipy netCDF4 pyyaml matplotlib h5py xarray dask[array] memory-profiler || rp_error pre_exec

# ------------------------------------------------------------------------------
# execute rank
$RP_PROF rank_start ""
mprof run -o mem_Grls_preprocess.dat "FittedISMIP_GrIS_preprocess.py" "--scenario" "ssp585" "--pipeline_id" "test.GrIS1f.FittedISMIP.GrIS" &

RP_EXEC_PID=$$
RP_RANK_PID=$!

wait $RP_RANK_PID
RP_RET=$?
$RP_PROF rank_stop "RP_EXEC_PID=$RP_EXEC_PID:RP_RANK_PID=$RP_RANK_PID"

# ------------------------------------------------------------------------------
# post-exec commands
$RP_PROF exec_post ""

# ------------------------------------------------------------------------------
$RP_PROF exec_stop ""
exit $RP_RET

# ------------------------------------------------------------------------------

