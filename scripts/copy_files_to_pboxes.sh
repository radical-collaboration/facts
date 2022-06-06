#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=copy_pb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --time=24:00:00
#SBATCH --output=out_copy_files_to_pboxes.txt
#SBATCH --error=error_copy_files_to_pboxes.txt
#SBATCH --export=ALL


# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base

echo "PB 1e"
bash copy_files_to_pb1e_pbox.sh
bash copy_files_to_pb1e_pbox_global.sh
bash copy_files_to_pb1e_pbox_rates.sh
bash copy_files_to_pb1e_pbox_rates_global.sh

echo "PB 2e"
bash copy_files_to_pb2e_pbox.sh
bash copy_files_to_pb2e_pbox_global.sh
bash copy_files_to_pb2e_pbox_rates.sh
bash copy_files_to_pb2e_pbox_rates_global.sh

echo "PB 1f"
bash copy_files_to_pb1f_pbox.sh
bash copy_files_to_pb1f_pbox_global.sh
bash copy_files_to_pb1f_pbox_rates.sh
bash copy_files_to_pb1f_pbox_rates_global.sh

echo "PB 2f"
bash copy_files_to_pb2f_pbox.sh
bash copy_files_to_pb2f_pbox_global.sh
bash copy_files_to_pb2f_pbox_rates.sh
bash copy_files_to_pb2f_pbox_rates_global.sh
