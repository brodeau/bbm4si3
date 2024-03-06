[![DOI](https://zenodo.org/badge/738921416.svg)](https://zenodo.org/badge/latestdoi/738921416)

# bbm4si3

A set of new and modified NEMO Fortran sources required to allow the sea-ice
component of NEMO, namely SI3, to use the new BBM brittle rheology
of [Òlason *et al.*, 2023](https://doi.org/10.1029/2021MS002685) in place of the
viscous-plastic rheologies available in SI3.

So far, these files are compatible with all `4.2.*` NEMO releases, visit
https://forge.nemo-ocean.eu/nemo to obtain the official NEMO code.


### How to test?

Simply copy all the `*.F90` and `*.h90` files found into `nemo_4.2/MY_SRC/` into
the `MY_SRC` directory of the NEMO configuration you wish to run (*i.e.* under
`<NEMO_REPO>/cfgs/<NEMO_CONFIG>/MY_SRC/`).  Then, compile NEMO as usual, the sources
found in the `MY_SRC` will overwrite those of the NEMO version you are using.


No-slip conditions must be used at the coast! Not only for SI3 but also for OPA !!!

L. Brodeau, January 2024
