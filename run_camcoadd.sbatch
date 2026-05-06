#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --job-name=all-64
#SBATCH --array=1
#SBATCH --output=sbatch_files/%x-%a.out
#SBATCH --error=sbatch_files/%x-%a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=dr.raga.pucha@gmail.com

./run_camcoadd.sh ${SLURM_ARRAY_TASK_ID} &
wait