#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=copy_wf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=2:00:00
#SBATCH --output=out_copy_files_to_workflows.txt
#SBATCH --error=error_copy_files_to_workflows.txt
#SBATCH --export=ALL


echo "WF 1e"
bash copy_files_to_wf1e_workflow.sh
bash copy_files_to_wf1e_workflow_global.sh

echo "WF 1f"
bash copy_files_to_wf1f_workflow.sh
bash copy_files_to_wf1f_workflow_global.sh

echo "WF 2e"
bash copy_files_to_wf2e_workflow.sh
bash copy_files_to_wf2e_workflow_global.sh

echo "WF 2f"
bash copy_files_to_wf2f_workflow.sh
bash copy_files_to_wf2f_workflow_global.sh

echo "WF 3e"
bash copy_files_to_wf3e_workflow.sh
bash copy_files_to_wf3e_workflow_global.sh

echo "WF 3f"
bash copy_files_to_wf3f_workflow.sh
bash copy_files_to_wf3f_workflow_global.sh

echo "WF 4"
bash copy_files_to_wf4_workflow.sh
bash copy_files_to_wf4_workflow_global.sh
