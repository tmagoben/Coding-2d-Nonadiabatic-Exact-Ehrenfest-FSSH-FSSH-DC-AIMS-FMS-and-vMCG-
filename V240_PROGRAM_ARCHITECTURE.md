# v0.24.0 program architecture

```text
EXTERNAL EVIDENCE PLANE                 CALLER TRUST PLANE

OpenMolcas 26.06 records                exact protocol fingerprint
  input / output / rassi.h5             exact SOC convention
  exported H_sf, H_SOC                  parser + exporter identity
  cross-geometry overlaps               manifest + environment SHA-256
             |                                      |
             +------------------+-------------------+
                                v
                    strict v0.24.0 bundle parser
                    - exact file inventory
                    - path/symlink rejection
                    - native/export digest binding
                    - fixture/external separation
                                |
               +----------------+----------------+
               |                                 |
               v                                 v
  transported derivative audit       independent validation audit
  9 coordinates x 3 steps             reference backend
  H_sf and H_SOC separated             basis/method ladders
  raw overlap -> polar transport       rigid frame invariance
  complete S/T manifold                state-tracking quality
               |                                 |
               +----------------+----------------+
                                v
                    fail-closed admission audit
                     /                      \
             fixture/protocol          admitted external
             diagnostics only          snapshot dynamics
                     |                       |
                     v                       v
             never production       frozen-nuclei unitary
                dynamics            propagation + restart
```

The existing Gaussian block dynamics, spin-free PySCF runtime, analytic SOC models,
replay v2, NAC convention, finite-manifold transport, and complete-multiplet controls
remain inherited below this new evidence plane.

