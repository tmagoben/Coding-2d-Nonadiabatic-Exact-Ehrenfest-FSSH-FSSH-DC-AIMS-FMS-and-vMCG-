# v0.25.1 release notes

Release date: 2026-08-24

v0.25.1 is the frozen-width multi-Gaussian TDVP metric-layer release. It replaces
the v0.25.0 single-canonical-packet restriction with multiple analytically coupled
Gaussian packets, while keeping the scientific domain narrow enough that every
matrix element and solver receipt can be independently recomputed.

## Added

- `FrozenGaussianSpinorStateV251`: packet-specific centers, momenta, frozen positive
  widths, and a complete electronic spinor on every packet.
- `QuadraticSpinHamiltonianV251`: a fixed-frame, one-dimensional Hermitian
  `H(x)=H0+x H1+x^2 H2` contract with positive constant nuclear mass.
- Exact analytic unequal-width Gaussian moments through degree three, sufficient for
  overlaps, kinetic energy, quadratic potentials, tangent metrics, and TDVP forcing.
- A real-parameter McLachlan system `G theta_dot=b` over coefficient real/imaginary
  parts and every packet center/momentum.
- Full-SVD metric receipts with absolute/relative rank cutoffs, retained condition
  number, projected null-space forcing, pseudoinverse residual, and velocity norm.
- A fully implicit midpoint nonlinear solve with an explicit predictor, solver
  status/counts, and a recomputable endpoint residual.
- Complete endpoint receipts binding model, settings, widths, midpoint parameters,
  metric/RHS/velocity/SVD, norm, energy, signed time, and nonlinear convergence.
- 55 deterministic numerical gates and 20 adversarial/core gates, yielding 75 new
  and 535 cumulative release gates.

## Validated reductions and covariances

- complete even singlet/triplet and odd Kramers-doublet analytic SOC spinors;
- exact continuous one-packet harmonic reductions `qdot=p/m` and `pdot=-kq`;
- compatible metric rank deficiency for exactly duplicate packets;
- signed-step implicit-midpoint reversal;
- arbitrary Gaussian relabeling;
- constant complex unitary electronic-frame transformations;
- exact zero-SOC enabled/disabled equivalence;
- second-order timestep refinement.

## Claim boundary

v0.25.1 is a one-dimensional, fixed-width, fixed-electronic-frame, quadratic-model
multi-Gaussian TDVP layer. It does not validate adaptive widths, packet spawning or
pruning, coordinate-dependent electronic frames in this solver, multidimensional
multi-Gaussian motion, real PySCF SOC trajectories, or general ab-initio dynamics
accuracy. It is therefore a foundation for later general variational dynamics, not
a production AIMS or adaptive vMCG implementation.

Recompute with:

```bash
python examples/133_recompute_v0251_multigaussian_tdvp.py
python examples/134_recompute_v0251_campaign.py
```

