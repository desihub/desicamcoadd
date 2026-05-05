# Script to coadd cameras for DESI-fuji
# A. Bolton, Jun 2022

# Imports and setups:
import os
import numpy as np
#import healpy as hp
from glob import glob
#import fitsio
from collections import defaultdict
#import desispec.io
import time
#import pandas as pd

#import matplotlib.pyplot as plt
#import h5py
#from astropy.table import Table
#from astropy.io import fits
#from redrock import templates
#from redrock import rebin
#from redrock import targets
#from redrock import zscan
#from desispec import coaddition
#from prospect import viewer

import camcoadd_funcs as ccf

#%set_env DESI_SPECTRO_REDUX=/global/cfs/cdirs/desi/spectro/redux
#%set_env SPECPROD=fuji
desi_spectro_redux='/global/cfs/cdirs/desi/spectro/redux'
specprod='fuji'

# Set the output directory:
outdir = '/global/cfs/cdirs/desi/science/gqp/camcoadd/' + specprod

# Construct the top-level healpix data directory:
#topdir = os.getenv("DESI_SPECTRO_REDUX") + "/" + os.getenv("SPECPROD") + '/healpix'
topdir = desi_spectro_redux + "/" + specprod + '/healpix'

# Specify the list of surveys to work on:
survey_list = ['cmx', 'special', 'sv1', 'sv2', 'sv3']
# Specify the list of programs to work on:
program_list = ['backup', 'bright', 'dark', 'other']
# Figure out which of these survey-plus-program combos exist:
s_p_combos = []
for this_survey in survey_list:
        for this_program in program_list:
                this_combo = this_survey + '/' + this_program
                if os.path.exists(topdir + '/' + this_combo):
                        s_p_combos.append(this_combo)

# Get a full list of lowest-level healpix directories:
healdir_list = []
for this_combo in s_p_combos:
        healdir_list += glob(topdir + '/' + this_combo + '/*/*')

# Translate into target output directories:
healout_list = [this_dir.replace(topdir, outdir + '/healpix', 1) for this_dir in healdir_list]
ndir = len(healdir_list)

print('Found ', ndir, ' input directories')

for i in range(ndir):
#for i in range(3):
        print(i, healdir_list[i])
        this_ofile = ccf.desi_camcoadd_healpix(healdir_list[i], healout_list[i])
