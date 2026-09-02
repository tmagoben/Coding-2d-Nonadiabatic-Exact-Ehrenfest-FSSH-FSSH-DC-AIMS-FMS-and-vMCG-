# v0.26.0 Validation

## Gate inventory

| Layer | Gates |
|---|---:|
| Inherited through v0.25.3 | 715 |
| New scientific validation | 80 |
| New adversarial/core controls | 30 |
| New in v0.26.0 | 110 |
| Cumulative | **825** |

## Principal numerical results

| Quantity | Result | Gate |
|---|---:|---:|
| Exact-grid norm drift | `0.0` | `<2e-12` |
| Exact-grid reversal error | `4.50e-15` | `<2e-11` |
| Edge-strip probability | `8.23e-25` | `<1e-12` |
| Observed grid time order | `2.07`, `2.32` | `>1.8` |
| Independent overlap quadrature error | `1.27e-14` | `<2e-10` |
| Independent Hamiltonian quadrature error | `3.47e-17` | `<2e-9` |
| 1D metric reduction error | `9.70e-17` | `<3e-13` |
| 1D velocity reduction error | `1.47e-16` | `<3e-13` |
| TDVP signed reversal error | `9.26e-12` | `<2e-8` |
| Constant-gauge velocity error | `1.82e-10` | `<2e-8` |
| Packet-permutation velocity error | `9.03e-11` | `<2e-8` |
| Signed-axis velocity error | `9.03e-11` | `<2e-8` |
| One-packet exact wavefunction error | `1.20e-3` | `<3e-3` |
| Controlled adaptive error | `6.73e-4` | `<3e-3` |
| Adaptive/one-packet error ratio | `0.561` | `<0.8` |

## Independent validations

- FFT-grid overlap and Hamiltonian matrix elements versus analytic moments.
- Exact-grid Strang timestep refinement against a finer independent trajectory.
- Exact one-dimensional reduction against v0.25.2 matrices, metric, RHS, and
  velocity.
- Exact-grid versus Gaussian wavefunctions on a common normalized grid.
- Reduced electronic density comparison independent of global wavefunction phase.

## Symmetry validations

- Hermiticity at arbitrary coordinates.
- Analytic potential derivatives versus centred finite differences.
- Coordinate transformation of the model potential.
- Complete Kramers time-reversal symmetry and pair degeneracy.
- Complete singlet and triplet projector ranks.
- Constant complex electronic gauge covariance.
- Packet permutation covariance.
- Signed coordinate-permutation covariance.
- Zero-SOC equivalence.

## Lifecycle validations

- all signed coordinate and momentum candidates;
- duplicate-candidate zero novelty and rank rejection;
- highest admitted residual score selection;
- exact zero-coefficient spawn projection;
- stable monotone IDs and zero newborn age;
- low-population aged pruning through projection gates;
- high-overlap merge-to-survivor through projection gates;
- dormant shape freezing with active coefficient transfer;
- population-plus-metric activation;
- no-event exact reduction to the base TDVP step;
- packet-cap fail-closed behavior.

## Negative/adversarial validations

Invalid masses, Hermiticity, projector resolution, grids, timesteps, widths,
coordinate transformations, active masks, projection rank, packet metadata, serial
identity, nonlinear convergence, and tampered receipts are rejected.

## Interpretation

The short-time adaptive improvement demonstrates that residual-driven basis growth
adds useful missing variational content on the frozen model.  It does not establish
long-time convergence, molecular accuracy, or equivalence to production AIMS/vMCG.
