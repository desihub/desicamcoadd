#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=04:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --job-name=desicamcoadd-64
#SBATCH --array=0-63
#SBATCH --output=/pscratch/sd/b/bweaver/%x-%j-%2a.out
#SBATCH --error=/pscratch/sd/b/bweaver/%x-%j-%2a.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=[${USER}@nersc.gov]
source /global/common/software/desi/desi_environment.sh main
# module swap desispec/0.99.99
export DESICAMCOADD=${HOME}/Documents/Code/git/desihub/desicamcoadd
export PYTHONPATH=${DESICAMCOADD}:${PYTHONPATH}
export SPECPROD=loa
cd ${DESICAMCOADD}
srun -N 1 -n 1 python desicamcoadd.py --batch ${SLURM_ARRAY_TASK_ID} --verbose ${SCRATCH}
