#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --job-name=od8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128000
#SBATCH --time=172:00:00
#SBATCH --output=out_od8.txt
#SBATCH --error=error_od8.txt
#SBATCH --export=ALL

cd /scratch/ggg46/localize_projections/oceandynamics

# Activate the conda environment
source /home/ggg46/.bashrc
conda activate base


LFILE="../location_lists/location_list.lst"

# Original seed - 1463
# Seed to correlate input temperatures - 4321


echo "OD - tlim4.0win0.25"
srun python tlm_preprocess_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --scenario tlim4.0win0.25 --baseyear 2005 --pyear_start 2020 --pyear_end 2300 --pyear_step 10 --model_dir ./data/cmip6 --locationfile ${LFILE}
srun python tlm_fit_oceandynamics.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25
srun python tlm_project_oceandynamics_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --nsamps 20000 --seed 4321 --tlm 1
srun python tlm_postprocess_oceandynamics_dev_dev.py --pipeline_id oceandynamics-tlm-oceandynamics-tlim4.0win0.25 --chunksize 662 --nsamps 20000 --seed 4321


mv *tlim4.0win0.25*globalsl* ../output_global
mv *tlim4.0win0.25*localsl* ../output_local

#rm *tlim4.0win0.25*.pkl
