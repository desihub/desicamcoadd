#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --job-name=all-64
#SBATCH --array=1-64
#SBATCH --output=sbatch_files/%x-%a.out
#SBATCH --error=sbatch_files/%x-%a.err

# Ensure subdirectory sbatch_files exist in the same directory as this file
# Otherwise, run mkdir -p sbatch_files in the terminal from this directory

chmod +x ./run_camcoadd.sh

./run_camcoadd.sh ${SLURM_ARRAY_TASK_ID} &
wait
