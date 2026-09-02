# v0.6 Build Validation Report

Validated on 2026-08-12.

## Compilation and regression suite

All repository Python source files compiled successfully.

```text
[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m.[0m[32m     [100%][0m
[32m[32m[1m68 passed[0m[32m in 2.25s[0m[0m
```

The cumulative suite includes all v0.1-v0.5 tests plus v0.6 tests for:

- global maximum-overlap root assignment;
- state sign/gauge correction;
- NAC transformation under root permutation and sign changes;
- ambiguity detection from weak overlaps and tied assignments;
- degenerate-subspace overlap singular values;
- unitary Procrustes/local-diabatic alignment;
- overlap-derived directional NACs;
- exact core+active CI determinant embedding;
- many-electron overlap with core-active cross-subspace mixing;
- root swaps/sign flips recovered from nonorthogonal many-electron overlaps;
- tracked backend state transformation and reset behavior;
- sequential tracked scans;
- order-independent interpolation of a completed 1D tracked scan.

## Synthetic crossing regression

```text
Raw overlap matrix:
[[ 0.05 -0.96]
 [ 0.94  0.03]]

Tracked -> raw permutation: [1 0]
Phase/sign corrections: [-1.+0.j  1.+0.j]
Assigned positive overlaps: [0.96-0.j 0.94+0.j]

Raw energy order: [0.1 0.9]
Tracked identity order: [0.9 0.1]

Raw NAC:
[[ 0.   0.3]
 [-0.3  0. ]]
Tracked/gauge-corrected NAC:
[[ 0.   0.3]
 [-0.3  0. ]]
```

## Local-diabatic / overlap-NAC regression

```text
Electronic overlap matrix:
[[ 0.96891242+0.j -0.24740396+0.j]
 [ 0.24740396+0.j  0.96891242+0.j]]

Current-basis Procrustes rotation Q:
[[ 0.96891242+0.j  0.24740396+0.j]
 [-0.24740396+0.j  0.96891242+0.j]]

Overlap after local diabatic alignment:
[[ 1.00000000e+00+0.j -1.99535592e-17+0.j]
 [ 6.93916277e-18+0.j  1.00000000e+00+0.j]]

Principal-overlap singular values:
[1. 1.]

Directional NAC recovered from a small overlap step:
[[ 0. +0.j -0.7+0.j]
 [ 0.7+0.j  0. +0.j]]
Expected off-diagonal magnitude: 0.7
```

## PySCF runtime status

PySCF is not installed in the build environment.

The real PySCF v0.6 path is therefore validated through:

1. the v0.5 explicit PySCF call-contract regression tests;
2. the current official PySCF MCSCF, FCI-overlap, and cross-AO-overlap API contracts;
3. a deterministic fake nonorthogonal FCI overlap implementation used to validate the
   core+active state-overlap algebra;
4. backend-independent tracking, phase, NAC, degeneracy, and scan tests.

For real PySCF validation after installation:

```bash
pip install -e ".[pyscf]"

python examples/17_pyscf_tracked_lih_scan.py
python examples/18_pyscf_overlap_nac_check.py
```

The second example compares the finite-step overlap-derived directional derivative
coupling with the analytic PySCF NAC projected along the same LiH displacement.

## Scientific interpretation

Passing these tests demonstrates internal consistency of the tracking/gauge
implementation.

It does not make individual adiabatic states uniquely defined at an exact degeneracy.
When tracking is ambiguous, the default backend behavior is to raise, and the
subspace/projector diagnostics should be used instead.
