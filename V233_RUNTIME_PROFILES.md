# v0.23.3 runtime profiles

v0.23.3 separates two questions that earlier provenance records could blur.

## Release locked

The `release_locked` profile asks whether execution matches the canonical
validation environment closely enough to claim release-byte identity. It binds
CPython 3.12.13, Linux x86-64/little-endian, NumPy 2.5.2, SciPy 1.18.0, h5py
3.16.0, PySCF 2.13.1, verified PySCF RECORD content, the Python executable hash,
approved memory-probe modes, and one-thread numerical environment variables.

## Scientifically compatible

The `scientifically_compatible` profile asks whether an environment is inside the
declared scientific support window: Python 3.10 through 3.13, NumPy 1.24 through
2.x, SciPy 1.10 through 1.x, h5py 3.8 through 3.x, and exactly PySCF 2.13.1.
It does not claim byte identity.

Reports include every check and the profile fingerprint. A runtime can be
scientifically compatible while failing the release lock; the campaign includes
this as a required control.
