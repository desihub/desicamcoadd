#!/bin/bash

cd
source .bashrc
cd /global/cfs/cdirs/desicollab/science/gqp/camcoadd/loa_code/

# Latest (not tagged)
# source /global/cfs/cdirs/desi/software/desi_environment.sh main
# Tagged for Loa
source $CFS/desi/software/desi_environment.sh 25.3

# can manually set template path to ensure correct templates get used
# Stephen Bailey confirmed Iron and Loa templates:
# DR1/Iron desimodules/23.1 redrock/0.17.0 redrock-templates/0.8 prospect/1.2.4
# DR2/Loa desimodules/24.11 redrock/0.20.3 redrock-templates/0.9 prospect/1.3.3
# export RR_TEMPLATE_DIR="/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/code/redrock-templates/0.9.1"

python run_camcoadd.py camcoadd_exec_$1.py
