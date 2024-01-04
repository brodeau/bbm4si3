# bbm4si3

Set of new and modified NEMO Fortran sources required to allow the sea-ice
component of NEMO, SI3, to use the new BBM brittle rheology
of [Òlason et al., 2023](https://doi.org/10.1029/2021MS002685) rather than the
default viscous-plastic rheologies of SI3.

So far, these files are compatible with all `4.2.*` NEMO releases, visit
https://forge.nemo-ocean.eu/nemo to obtain the official NEMO code.


### How to test?

Simply copy all the `*.F90` and `*.h90` files found into `nemo_4.2/MY_SRC/` into
the `MY_SRC` directory of the NEMO configuration you wish to run (*i.e.* under
`<NEMO_REPO>/cfgs/<NEMO_CONFIG>/MY_SRC/`).  Then compile as usual, the sources
found in the `MY_SRC` will overwrite those of the NEMO version you use.

