#!/bin/bash

export RP_RESOURCE_SANDBOX="/opt"

## declare array of workflows from worflows.yml
# declare -a arr=("wf1e" "wf1f" "wf2e" "wf2f" "wf3e" "wf3f" "wf4")


# for workflow in "wf1e" "wf1f" "wf2e" "wf2f" "wf3e" "wf3f" "wf4" ;             # problem  "wf3e"
for workflow in "wf1e" "wf1f" "wf2e" "wf2f" "wf3f" "wf4" ;
    do
        mprof run -o mem_total_$workflow.dat "total_workflow.py" "--directory" "$RP_RESOURCE_SANDBOX/to_total" "--workflows" "workflows.yml" "--workflow" $workflow "--scale" "local" "--experiment_name" "configTest" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10"
        # python3 total_workflow.py "--directory" "$RP_RESOURCE_SANDBOX/to_total" "--workflows" "workflows.yml" "--workflow" $workflow "--scale" "global" "--experiment_name" "configTest" "--pyear_start" "2020" "--pyear_end" "2150" "--pyear_step" "10"
    done

echo "=== FINISHED TOTALING STEP ==="
