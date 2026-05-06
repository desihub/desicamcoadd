"""
This is a wrapper function to run the `desi_camcoadd_healpix` function with multiprocessing.

Author: Ragadeepika Pucha, Adam Bolton, Stéphanie Juneau
Version: 2025 April 18
"""

####################################################################################################
import sys
from camcoadd_funcs import desi_camcoadd_healpix
from multiprocessing import Pool

import time
####################################################################################################

## Input Filename
filename = str(sys.argv[1])

f = open(filename, 'r')
lines = f.readlines()[1:]

in_args = [thisline.split("'")[1] for thisline in lines]
out_args = [thisline.split("'")[3] for thisline in lines]

## Starting timer
start = time.time()

## Number of files
n_files = len(lines)
## Number of cores
n_cores = 64

args = [(in_args[idx], out_args[idx]) for idx in range(n_files)]
pool = Pool(processes = n_cores)
pool.starmap(desi_camcoadd_healpix, args)
pool.close()
pool.join()

end = time.time()
print (f'Time Taken for {n_files:4d} files with {n_cores:3d} cores: {round(end-start, 2)} seconds')

####################################################################################################