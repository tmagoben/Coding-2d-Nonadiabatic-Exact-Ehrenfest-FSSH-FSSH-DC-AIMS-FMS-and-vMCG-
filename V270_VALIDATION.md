# v0.27.0 Validation

## Scientific evidence: 100 gates

The new evidence is divided into six auditable groups:

| Group | Gates | Primary oracle |
|---|---:|---|
| Schemas and honest claim boundaries | 8 | Frozen constants and false claims |
| Symmetric/log-Euclidean matrix algebra | 12 | Isometries and centered finite differences |
| Correlated moments and matrix elements | 15 | Independent 192×192 FFT quadrature |
| v0.26.0/one-dimensional reductions | 11 | Inherited analytic implementation |
| Riccati dynamics and covariance | 35 | Exact matrix Riccati ODE, rotations, reflections, gauge, permutations |
| Correlated basis lifecycle | 19 | Recomputed novelty, residual, projection, activation receipts |

The total is 100; the table combines Riccati dynamics and covariance as one
method block.

## Independent numerical results

- Maximum fourth-order moment error against direct grid quadrature:
  approximately `2.1e-16`.
- Maximum overlap-matrix grid error: approximately `1.6e-15`.
- Maximum Hamiltonian-matrix grid/FFT error: approximately `5.1e-17`.
- Riccati endpoint errors at `dt = 0.04, 0.02, 0.01`:
  approximately `7.60e-6`, `1.90e-6`, `4.75e-7`.
- Observed midpoint orders: approximately `1.9998` and `1.99995`.
- Arbitrary-rotation metric covariance error: below `4e-15`.
- Arbitrary-rotation midpoint endpoint error: below `7e-16`.
- Dormant full-shape drift after exact restoration: exactly zero.

## Adversarial controls: 35 gates

Controls reject disabled frozen contracts, live molecular SOC, invalid width domains,
non-SPD/nonsymmetric matrices, wrong projection shapes, duplicate projection bases,
tampered spawn-direction policies, invalid lifecycle metadata, forced nonlinear
nonconvergence, and modified metric/midpoint/projection/event receipts.

## Cumulative campaign

- inherited v0.26.0: 825 gates;
- v0.27.0 scientific evidence: 100 gates;
- v0.27.0 adversarial/core controls: 35 gates;
- total: **960/960**.

The campaign does not convert analytic benchmark success into a molecular accuracy
claim. Live PySCF SOC trajectories and coordinate-dependent electronic frames remain
closed.
