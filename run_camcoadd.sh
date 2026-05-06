#!/bin/bash

cd
source .bashrc
cd /global/cfs/cdirs/desicollab/science/gqp/camcoadd/loa_code/

# Latest (not tagged)
# source /global/cfs/cdirs/desi/software/desi_environment.sh main
# Tagged for Loa
source $CFS/desi/software/desi_environment.sh 25.3

# manually set template path to ensure correct templates get used
export RR_TEMPLATE_DIR="/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/code/redrock-templates/main"

python run_camcoadd.py camcoadd_exec_$1.py