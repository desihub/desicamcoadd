# Licensed under a 3-clause BSD style license - see LICENSE.md.
"""
desicamcoadd
============

Functions for coadding DESI spectra across cameras.

-- A. Bolton, Jun 2022
-- Substantially revised, B. A. Weaver, Nov 2026
"""
import os
import sys
import time
import json
import argparse
from glob import glob
import numpy as np
from astropy.io import fits
# import healpy as hp
# import fitsio
# from collections import defaultdict
from desispec import coaddition
from desiutil.log import get_logger, DEBUG
import desispec.io
from prospect import viewer


def compute_resolution_sigma(spec_struc):
    """Compute wavelength-dependent "sigma" of line-spread function, in units of the wavelength baseline.

    **CURRENTLY ONLY WORKS FOR DATA COADDED ACROSS CAMERAS!!** (with 'brz' as the only band)

    Parameters
    ----------
    spec_struc : :class:`desispec.spectra.Spectra`
        A ``Spectra`` container object produced by :func:`~desispec.coaddition.coadd_cameras`.

    Returns
    -------
    :class:`numpy.ndarray`
        An array of wavelength dispersion for each spectral element in `spec_struc`.
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
    rsigma *= dwave  # using numpy broadcasting
    return rsigma


def desi_camcoadd_healpix(survey, program, healpix, output_root):
    """Function to coadd DESI healpix-format spectra across cameras.

    Parameters
    ----------
    survey : :class:`str`
        The survey name.
    program : :class:`str`
        The program name.
    healpix : :class:`str` or :class:`int`
        The healpix number.
    output_root : :class:`str`
        Write output file(s) to this directory.

    Returns
    -------
    :class:`str`
        Filename (with dirpath) of the written-out file, or ``None`` if one
        or more input files were missing.
    """
    log = get_logger()
    #
    # Construct filenames for input and output files
    #
    cfile, cfile_ok = desispec.io.findfile('coadd', survey=survey, faprogram=program, healpix=healpix, readonly=True, return_exists=True)
    zfile, zfile_ok = desispec.io.findfile('redrock', survey=survey, faprogram=program, healpix=healpix, readonly=True, return_exists=True)
    mfile, mfile_ok = desispec.io.findfile('rrmodel', survey=survey, faprogram=program, healpix=healpix, readonly=True, return_exists=True)
    ofile = os.path.join(output_root, survey, program, str(int(healpix) // 100), str(healpix), f"camcoadd-{survey}-{program}-{healpix}.fits")
    if (cfile_ok and (mfile_ok or zfile_ok)):
        # Adding time calculation - RP
        start_time = time.time()
        # Get spectrum & redshift structures:
        spec_struc = desispec.io.read_spectra(cfile)

        # Adding time calculation - RP
        step1_time = time.time()
        log.info('Time for reading files: %d seconds.', round(step1_time - start_time))

        # Coadd across cameras:
        spec_struc_2 = coaddition.coadd_cameras(spec_struc)
        # Compute approximating LSF sigma:
        rsigma = compute_resolution_sigma(spec_struc_2)
        # Generate model:
        if mfile_ok:
            # TO DO
            model_wave = dict()
            model_flux = dict()
            with fits.open(mfile) as hdulist:
                for band in 'BRZ':
                    model_wave[band.lower()] = hdulist[f"{band}_WAVELENGTH"].data
                    model_flux[band.lower()] = hdulist[f"{band}_MODEL"].data
            #
            # This is adapted from prospect.viewer.create_model.
            #
            keep = dict()
            meanwaves = [np.mean(model_wave[band]) for band in 'brz']  # sort bands by increasing wavelength
            sorted_bands = [x for _,x in sorted(zip(meanwaves, 'brz'))]
            for i_band, band in enumerate(sorted_bands):
                wavecut_low = 0
                wavecut_up = 1.e10
                if i_band > 0:
                    band_low = sorted_bands[i_band - 1]
                    wavecut_low = 0.5*(model_wave[band_low][-1] + model_wave[band][0])
                if i_band < len(sorted_bands) - 1:
                    band_up = sorted_bands[i_band + 1]
                    wavecut_up = 0.5*(model_wave[band][-1] + model_wave[band_up][0])
                keep[band] = (model_wave[band] > wavecut_low) & (model_wave[band] < wavecut_up)
            wavemodel = np.concatenate([model_wave[band][keep[band]] for band in sorted_bands])
            fluxmodel = np.concatenate([model_flux[band][:, keep[band]] for band in sorted_bands],
                                       axis=1)
        else:
            zstruc = fits.getdata(zfile, 1)
            wavemodel, fluxmodel = viewer.create_model(spec_struc_2, zstruc)

        # Append model (& eventually other stuff) to new structure using "extra" field:
        spec_struc_2.extra = {spec_struc_2.bands[0]: {'model': fluxmodel.copy(),
                                                      'sky': 0.*fluxmodel, # placeholder for sky
                                                      'wavedisp': rsigma}}

        # Adding time calculation - RP
        step2_time = time.time()
        log.info('Time for Coadding Spectra: %d seconds.', round(step2_time - step1_time))

        # Write out to file:
        log.debug("os.makedirs('%d', exist_ok=True)", os.path.dirname(ofile))
        os.makedirs(os.path.dirname(ofile), exist_ok=True)
        desispec.io.write_spectra(ofile, spec_struc_2)

        # Adding time calculation - RP
        end_time = time.time()
        log.info('Time for Writing Spectra: %d seconds,', round(end_time - step2_time))

        return ofile
    else:
        log.error('Missing one or more input files for %s, %s, %s.', survey, program, healpix)
        return None


def initialize_input_healpix(batch_file, specprod=None, flavour='healpix'):
    """Create a list of input healpix coadds, and organize them into batches.

    Parameters
    ----------
    batch_file : :class:`str`
        Write the input batches to this file.
    specprod : :class:`str`, optional
        Name of the specprod. If not defined, use :envvar:`SPECPROD`.
    flavour : :class:`str`, optional
        Type of coadd to search for: ``healpix`` for ``iron``, ``loa``, but
        may change in the future.

    Returns
    -------
    :class:`int`
        The number of input directories found.
    """
    log = get_logger()
    if specprod is None:
        specprod = os.environ['SPECPROD']
    search_dir = os.path.join(os.environ['DESI_SPECTRO_REDUX'], specprod, flavour)
    s_p_combos = []
    for survey in ('cmx', 'special', 'sv1', 'sv2', 'sv3', 'main'):
        for program in ('backup', 'bright', 'dark', 'other'):
            if os.path.exists(os.path.join(search_dir, survey, program)):
                s_p_combos.append((survey, program))
    #
    # Find all directories.
    #
    healpix_directories = []
    for combo in s_p_combos:
        log.debug("glob(os.path.join('%s', '%s', '%s' ,'*' , '*'))", search_dir, combo[0], combo[1])
        healpix_directories += glob(os.path.join(search_dir, combo[0], combo[1], '*', '*'))
    #
    # Convert to survey, program, healpix.
    #
    sph = [tuple(d.rsplit('/', 4)[1:]) for d in healpix_directories]
    healpix_sizes = []
    for survey, program, hpxgroup, healpix in sph:
        cfile, cfile_ok = desispec.io.findfile('coadd', survey=survey,
                                               faprogram=program, healpix=healpix,
                                               readonly=True, return_exists=True)
        if cfile_ok:
            healpix_sizes.append(os.path.getsize(cfile))
        else:
            healpix_sizes.append(0)
    #
    # Find the number of objects in each coadd.
    # Then batch up files so that the number of total objects, rather than
    # the number of files is roughly equal in each batch.
    # As a first pass though, try just batching by total size.
    #
    with open(batch_file, 'w') as JSON:
        json.dump(sph, JSON)
    return len(healpix_directories)


def _options():
    """Parse command-line options.
    """
    prsr = argparse.ArgumentParser()
    prsr.add_argument('-b', '--batch', action='store', metavar='N', type=int,
                      help='Process batch N')
    prsr.add_argument('-i', '--initialize', action='store_true',
                      help='Find coadd files and split them into batches.')
    prsr.add_argument('-v', '--verbose', action='store_true',
                      help='Print extra debugging information.')
    prsr.add_argument('output', metavar='DIR', help='Write output to DIR.')
    return prsr.parse_args()


def main():
    """Entry-point for command-line scripts.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    options = _options()
    if options.verbose:
        log = get_logger(DEBUG)
    else:
        log = get_logger()
    specprod = os.environ['SPECPROD']
    if options.initialize:
        n = initialize_input_healpix(os.path.join(options.output, f'desicamcoadd_{specprod}.json'))
        log.info("Found %d coadd files.", n)
    else:
        if options.batch is None:
            log.critical('Batch number not specified!')
            return 1
        else:
            log.warning('Processing mode is a work in progress.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
