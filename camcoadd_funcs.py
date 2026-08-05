# camcoadd_funcs.py
# functions for coadding DESI spectra across cameras
# A. Bolton, Jun 2022

# Imports and setups:
import os
import numpy as np
import healpy as hp
from glob import glob
import fitsio
from collections import defaultdict
import desispec.io
import time
import pandas as pd

import matplotlib.pyplot as plt
import h5py
from astropy.table import Table
from astropy.io import fits
from redrock import templates
from redrock import rebin
from redrock import targets
from redrock import zscan
from desispec import coaddition
from prospect import viewer

def compute_resolution_sigma(spec_struc):
    """Function to compute wavelength-dependent "sigma" of line-spread function,
    in units of the wavelength baseline
    !!CURRENTLY ONLY WORKS FOR DATA COADDED ACROSS CAMERAS!! (with 'brz' as the only band)
    """
    # Get dimensionality of resolution data array:
    nspec, nband, npix = spec_struc.resolution_data['brz'].shape
    # Create arrays for calculating LSF moments:
    xband = np.arange(float(nband))
    xband -= xband.mean()
    xfull = np.outer(xband, np.ones(npix))
    rsigma = np.full((nspec, npix), 0.)
    # Loop over spectra to compute dispersion values, initially in units of pixels:
    for ispec in range(nspec):
        rnorm = spec_struc.resolution_data['brz'][ispec].sum(0)
        rmask = rnorm > 0
        rnorm_inv = np.full(npix, 0.)
        rnorm_inv[rmask] = 1. / rnorm[rmask]
        xres = rnorm_inv * (spec_struc.resolution_data['brz'][ispec] * xfull).sum(0)
        x2res = rnorm_inv * (spec_struc.resolution_data['brz'][ispec] * xfull**2).sum(0)
        rsigma[ispec] = np.sqrt(np.abs(x2res - xres**2))
    # Convert from pixels to delta-wavelength:
    dwave = 0. * spec_struc.wave['brz']
    dwave[1:-1] = 0.5 * (spec_struc.wave['brz'][2:] - spec_struc.wave['brz'][:-2])
    dwave[0] = dwave[1]
    dwave[-1] = dwave[-2]
    rsigma *= dwave # using numpy broadcasting
    return rsigma

def desi_camcoadd_healpix(inputdir, outputdir):
    """Function to coadd DESI healpix-format spectra across cameras.
    Finds input files in inputdir and writes output files to outputdir.
    Returns filename (with dirpath) of the written-out file.
    """
    # Parse the inputdir into individual survey, program, and healpix strings:
    hp_string = inputdir.split('/')[-1] # healpix string
    hs_string = inputdir.split('/')[-2] # healpix substring
    pr_string = inputdir.split('/')[-3] # program string
    su_string = inputdir.split('/')[-4] # survey string
    # Construct filenames for input and output files:
    fn_string = su_string + "-" + pr_string + "-" + hp_string
    cfile = inputdir + "/" + "coadd-" + fn_string + ".fits" # input coadd file
    zfile = inputdir + "/" + "redrock-" + fn_string + ".fits" # input file for redshift mode data
    ofile = outputdir + "/" + "camcoadd-" + fn_string + '.fits'
    # sfile = inputdir + "/" + "spectra-" + fn_string + ".fits" # may need in future
    # rfile = inputdir + "/" + "rrdetails-" + fn_string + ".h5" # may need in future
    # Test for existence of input files:
    if (os.path.exists(cfile) and os.path.exists(zfile)):
        ## Adding time calculation - RP
        start_time = time.time()
        # Get spectrum & redshift structures:
        spec_struc = desispec.io.read_spectra(cfile)
        zstruc = fits.getdata(zfile,1)

        ## Adding time calculation - RP
        step1_time = time.time()
        print ('Time for reading files: ', round(step1_time - start_time), 'seconds')
        
        # Coadd across cameras:
        spec_struc_2 = coaddition.coadd_cameras(spec_struc)
        # Compute approximating LSF sigma:
        rsigma = compute_resolution_sigma(spec_struc_2)
        # Generate model:
        wavemodel, fluxmodel = viewer.create_model(spec_struc_2, zstruc)
        # Append model (& eventually other stuff) to new structure using "extra" field:
        spec_struc_2.extra = {spec_struc_2.bands[0]: {'model': fluxmodel.copy(),
                                                      'sky': 0.*fluxmodel, # placeholder for sky
                                                      'wavedisp': rsigma}}

        ## Adding time calculation - RP
        step2_time = time.time()
        print ('Time for Coadding Spectra: ', round(step2_time - step1_time), 'seconds')
        
        # Write out to file:
        os.makedirs(outputdir, exist_ok=True)
        desispec.io.write_spectra(ofile, spec_struc_2)

        ## Adding time calculation - RP
        end_time = time.time()
        print ('Time for Writing Spectra: ', round(end_time - step2_time), 'seconds')
        
        return ofile
    else:
        print('Missing an input file')
        return 0
