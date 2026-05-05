# desicamcoadd

Scripts to automate coaddition of DESI spectra across cameras.

## Operation

This script will loop over all healpix in a specprod (*e.g.* `fuji`)
and produce `camcoadd` files.

```shell
$ python camcoadd_fuji.py
```

## Warnings

* Ensure that the script is run in a Python environment consistent with the
  environment used to produce the specprod.
* The coaddition of sky spectra is not yet supported. The `camcoadd` files
  typically will contain a `BRZ_SKY` HDU, but all values will be set to zero.

## License

desicamcoadd is free software licensed under a 3-clause BSD-style license.
For details see the `LICENSE.md` file.

## Change Log

### 0.1.1 (skymodel test; c. 2022-09-19; tagged 2026-05-05)

This version represents an early attempt to add coadded sky spectra. This
appears to have resulted in `${DESI_ROOT}/science/gqp/camcoadd/test`, but
the (one) `camcoadd` file has zeroes for the `BRZ_SKY` HDU. Apparently,
the test was not successful.

### 0.1.0 (fuji; DESI EDR; c. 2022-06-07; tagged 2026-05-05)

Version used to generate the `camcoadd` files for `fuji`.
