# v0.25.0 release notes

Release date: 2026-08-24

v0.25.0 is the restricted time-dependent-variational SOC dynamics release. It
connects the already validated representation-neutral `H`, physical operator
derivative `K`, mass, and cross-geometry overlap contracts to a symmetric
single-packet propagator. It deliberately stops short of claiming the full coupled
multi-Gaussian TDVP.

## Added

- `CanonicalVariationalSOCStateV250` for constant-mass canonical nuclear coordinates
  and a normalized complete electronic spinor.
- A self-adjoint kick--drift--endpoint-Strang--kick step. Nuclear coordinates use
  velocity Verlet; electronic amplitudes use two endpoint Hamiltonian half steps and
  exact finite-manifold frame transport.
- SVD-based polar transport: for `O = U Sigma V^dagger`, amplitudes are transported
  with `W^dagger`, where `W = U V^dagger`; the raw singular values remain evidence.
- Step receipts that bind the endpoints to both mass matrices, both forces, both
  Hamiltonians, both physical derivative tensors, the raw overlap, its polar factor,
  quality policy, singular values, energies, and signed time step.
- Complete singlet/triplet and Kramers-doublet SOC trajectory tests, coordinate-
  dependent complex-gauge covariance, exact signed-step reversal, zero-SOC
  equivalence, norm preservation, and second-order timestep/energy behavior.
- Fail-closed controls for full multi-Gaussian TDVP, adaptive widths,
  coordinate-dependent mass under Verlet, static-only SOC snapshots, spectral
  expansion, manifold loss, non-SVD transport, and receipt tampering.
- 45 numerical validation gates and 15 core/adversarial gates, giving 60 new and 460
  cumulative release gates.

## Numerical decision

Velocity Verlet is the right restricted nuclear update here because the released
nuclear variables are canonical and their generalized mass is constant. It is not
declared suitable for the future full TDVP metric system. That system generally has
coupled, complex, noncanonical coordinates and should use an implicit midpoint or
discrete variational solve.

Polar decomposition and SVD are not competing choices. The polar factor is the
physical unitary transport; SVD is the robust algorithm used to compute it and to
expose retained-manifold singular values, condition number, and principal angles.

## Claim boundary

Validated in v0.25.0:

- the single canonical nuclear-packet / complete-spinor TDVP restriction;
- symmetric constant-mass Verlet plus endpoint-Strang coupling;
- SVD-computed unitary polar transport with raw-overlap quality gates;
- even- and odd-electron analytic SOC model trajectories;
- coordinate-dependent complex-gauge covariance and signed-step reversibility;
- second-order timestep and energy-error plateaus.

Not validated or admitted:

- full coupled multi-Gaussian or adaptive-width TDVP;
- plain Verlet on a general noncanonical variational manifold;
- coordinate-dependent generalized mass in the Verlet path;
- a trajectory-ready PySCF molecular-SOC provider;
- full Cartesian/analytic molecular SOC derivatives or a continuous physical
  molecular derivative connection;
- general ab-initio SOC dynamics accuracy.

Recompute with:

```bash
python examples/131_recompute_v0250_variational_soc.py
python examples/132_recompute_v0250_campaign.py
```
