#!/bin/bash

cd
source .bashrc
cd /global/cfs/cdirs/desicollab/science/gqp/camcoadd/loa_code/

# Latest (not tagged)
#source /global/cfs/cdirs/desi/software/desi_environment.sh main
# Tagged for Iron and Loa
source $CFS/desi/software/desi_environment.sh 25.3

python run_camcoadd.py camcoadd_exec_$1.py