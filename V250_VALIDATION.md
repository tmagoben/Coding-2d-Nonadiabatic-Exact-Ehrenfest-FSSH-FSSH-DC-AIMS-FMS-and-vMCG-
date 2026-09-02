# v0.25.0 validation

## Deterministic models

The numerical campaign uses both analytic SOC families introduced in v0.22.0:

- an even-electron singlet plus complete triplet;
- an odd-electron pair of complete Kramers doublets.

The odd model is also propagated in a coordinate-dependent complex unitary gauge.
A controlled `0.97 I` overlap separates the raw finite-manifold contraction from its
identity polar factor. A zero-SOC pair verifies exact enabled/disabled reduction.

## Numerical observations

- Even/odd maximum electronic-norm drift: below `2.2e-14`.
- Even/odd maximum absolute energy drift at `h=0.4` au: below `4.0e-10` hartree.
- Signed forward/reverse nuclear residuals: below `8.0e-15`.
- Signed forward/reverse spinor residuals: below `1.7e-14`.
- Coordinate-dependent complex-gauge spinor residual: below `1.7e-14`.
- Maximum unitary-polar residual: below `1.9e-15`.
- Minimum deliberately retained singular value: `0.97`.

At fixed final time `20` au, the step ladder `(0.8, 0.4, 0.2, 0.1)` gives successive
phase-aligned endpoint-change ratios near `0.25`; the energy-drift ratios also remain
near `0.25`. This is the expected second-order plateau.

## Acceptance structure

- 400 gates inherited unchanged from v0.24.2.
- 45 deterministic variational-SOC validation gates.
- 15 scope/adversarial/core gates.
- 60 new gates and 460 cumulative gates.

Negative controls reject requests for full TDVP, adaptive widths,
coordinate-dependent-mass Verlet, a non-SVD transport algorithm, spectrally
expansive/rank-lost overlaps, static-only SOC providers, and tampering with momentum,
spinor, polar transport, singular values, overlap metrics, or endpoint mass.

## Interpretation

This is implementation and structure validation on analytic SOC models. It is not a
general molecular-accuracy benchmark and does not turn the v0.24.1/v0.24.2 PySCF SOC
evidence into an arbitrary-geometry trajectory backend.
