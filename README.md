# desicamcoadd

Scripts to automate coaddition of DESI spectra across cameras.

## Operation

### 0. Environment

For new versions of `iron` and for `loa`, this is recommended:

```shell
source $CFS/desi/software/desi_environment.sh 25.3
```

### 1. Preliminary tests

`loa_dev_test_01.ipynb`: Preliminary tests from `iron` file copied over and modified for `loa`.

### 2. Cores and files per batch

Conduct tests on the debug queue to optimize the number of cores & files per batch:

* Aim: take up to 20 min in debug queue (then request 1.5 times longer with real sbatch job: 30min for 20min in debug queue).
* Results: 128 cores and 96 cores was too large; ran into OOM issues. **64 cores is optimal**
  - 572 dirs with 64 cores takes 300sec (5min) so we can aim for ~2000 dirs per batch
  - Loa: 56780 directories -->  32 batches of ~1775 (or 28 batches of ~2028)

### 3. Set up batches

`dr2_camcoadd_setup.ipynb`: Notebook to examine the files, compute disk space, and group
lists of healpix directories in batches.

* Results: 56780 directories occupying 28.6 TB
& Goal: create 32 `camcoadd_exec_*.py` files (possibly after removing the ~13 files at >4 GB to avoid OOM issues)

### 4. Prepare for batch submission

When happy that all `camcoadd_exec_*.py` files exist and the `camcoadd_funcs.py` file is ready, check on:

* `run_camcoadd.py`: wrapper script; make sure it's set to use 64 cores
* `run_camcoadd.sh`: make sure it's starting from correct specprod_code (loa_code for DR2)
* `run_camcoadd.sbatch`:
  - edit the jobname if needed;
  - edit the queue between "debug" (max 5 processes) or "regular" (can launch all):
    `#SBATCH --qos=debug`  vs. `#SBATCH --qos=regular`
  - change `#SBATCH --array=1` to include the range of `camcoadd_exec_*.py` file numbers,
    *e.g.*, `#SBATCH --array=1-32` to run all 32 or could split between SJ (1-16) and AB (17-32)?
  - change the email address to receive notification of: job start, end, fail (if applicable)

### 5. Submit jobs

When everything's ready, here are commands to launch job, check queue, or cancel job:

```
sbatch run_camcoadd.sbatch  # to launch (edit filename if we split between _sj and _ab versions)
squeue --me     # for a person to only check their own jobs in the queue
scancel JOB-ID  # to cancel the job
```

### 6. QA checks.

After it's finished running, revisit later steps of `dr2_camcoadd_setup.ipynb` to perform some checks.


## Warnings

* Some bugs were fixed in `coadd_cameras` in Fall 2024. Some of those were present for `iron`.
  See [desihub/desispec#2377](https://github.com/desihub/desispec/pull/2377).
* Ensure that the script is run in a Python environment consistent with the
  environment used to produce the specprod.
* The coaddition of sky spectra is not yet supported. The `camcoadd` files
  typically will contain a `BRZ_SKY` HDU, but all values will be set to zero.

## License

desicamcoadd is free software licensed under a 3-clause BSD-style license.
For details see the `LICENSE.md` file.

## Change Log

### 2.0.0 (loa; DESI DR2; Work In Progress)

This will be the version used to generate the `camcoadd` files for `loa`.

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
