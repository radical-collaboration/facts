#!/bin/bash

# ---------------------------------------------------
# Usage
#   source submit_docker_experiment.sh <experiment_name>
#
# Example
#   source submit_docker_experiment.sh coupling.ssp585
# ---------------------------------------------------

# Experiment name passed as first argument
EXPERIMENT="${1:-}"

# Check that an experiment name was provided
if [ -z "$EXPERIMENT" ]; then
    echo "Usage: source submit_docker_experiment.sh <experiment_name>"
    return 1 2>/dev/null || exit 1
fi

#---------------------------------------------------
# Change the timezone
TZVAR="America/New_York"
#---------------------------------------------------
EXPPATH="experiments_ssisls/$EXPERIMENT"
LOGFILE="${EXPPATH}/time_${EXPERIMENT}_$(whoami).log"
#---------------------------------------------------
# Start timestamp
{
    echo "===== START: $(TZ="$TZVAR" date '+%Y-%m-%d %H:%M:%S') EST ====="
    
    START_TIME=$(TZ="$TZVAR" date +%s)

    # Run FACTS
    python runFACTS.py "$EXPPATH"

    END_TIME=$(TZ="$TZVAR" date +%s)
    DURATION=$((END_TIME - START_TIME))

    # Format duration as HH:MM:SS
    HOURS=$((DURATION / 3600))
    MINUTES=$(((DURATION % 3600) / 60))
    SECONDS=$((DURATION % 60))
    FORMATTED_DURATION=$(printf "%02d:%02d:%02d" $HOURS $MINUTES $SECONDS)

    echo "===== END: $(TZ="$TZVAR" date '+%Y-%m-%d %H:%M:%S') ====="
    echo "DURATION: $FORMATTED_DURATION (HH:MM:SS)"

} 2>&1 | tee "$LOGFILE"