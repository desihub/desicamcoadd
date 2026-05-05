# desicamcoadd

Scripts to automate coaddition of DESI spectra across cameras.

## Operation

The notebook `dr1_camcoadd_setup.ipynb` loops over all healpix in a specprod
(*e.g.* `iron`) and sets up a series of batch jobs to create the `camcoadd` files.

The batch jobs were executed on a perlmutter exclusive node, presumably in
interactive mode. Additional details are in the notebook.

After completion of the batch jobs, the notebook can be used for quality
assurance.

## Warnings

* Ensure that the script is run in a Python environment consistent with the
  environment used to produce the specprod.
* The coaddition of sky spectra is not yet supported. The `camcoadd` files
  typically will contain a `BRZ_SKY` HDU, but all values will be set to zero.

## License

desicamcoadd is free software licensed under a 3-clause BSD-style license.
For details see the `LICENSE.md` file.

## Change Log

### 1.0.0 (iron; DESI DR1; c. 2025-04-18; tagged 2026-05-05)

This is the version used to generate the `camcoadd` files for `iron`.

For `iron` the notebook `dr1_camcoadd_setup.ipynb` replaces `camcoadd_fuji.py`.
The file `camcoadd_funcs.py` was reverted to the original `fuji`/`0.1.0` version,
*i.e.* without the skymodel test modifications.

### 0.1.1 (skymodel test; c. 2022-09-19; tagged 2026-05-05)

This version represents an early attempt to add coadded sky spectra. This
appears to have resulted in `${DESI_ROOT}/science/gqp/camcoadd/test`, but
the (one) `camcoadd` file has zeroes for the `BRZ_SKY` HDU. Apparently,
the test was not successful.

### 0.1.0 (fuji; DESI EDR; c. 2022-06-07; tagged 2026-05-05)

This is the version used to generate the `camcoadd` files for `fuji`.
