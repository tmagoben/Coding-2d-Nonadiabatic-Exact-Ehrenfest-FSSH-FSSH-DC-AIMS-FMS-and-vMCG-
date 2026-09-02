# v0.13 Validation Contract

v0.13 is accepted only if the new residual-driven machinery is independently verified
and every retained v0.1-v0.12 regression continues to pass.

## 1. Initial projection residual

The initial residual is

$$
r_N=\Psi-P_N\Psi.
$$

For pure residual-greedy selection the test suite requires:

- every accepted Gaussian lowers the squared projection residual;
- the residual decreases monotonically with basis size;
- the analytically predicted one-candidate gain agrees with the actual reprojection
  reduction.

The acceptance check does not allow a candidate to be called useful merely because it
was selected by the search routine.

## 2. Candidate redundancy

For candidate overlap vector

$$
s_i=\langle g_i|g_c\rangle,
$$

the new independent nuclear norm is

$$
n_c=1-s^\dagger S^{-1}s.
$$

Candidates are rejected if this quantity is below the configured orthogonal-norm
floor.

The proposed expanded overlap condition number is also checked before accepting a
candidate.

## 3. Prepared/vectorized implementation

The release benchmark uses a precomputed Gaussian dictionary for speed.

A regression test compares the prepared/vectorized greedy builder with the slower
direct implementation on the same small dictionary.

They must choose the same residual-greedy candidates and reproduce the same final
projection residual within numerical tolerance.

## 4. Observable-aware screen

The density-aware v0.13 selector is not allowed to search the full dictionary directly
by a reduced observable.

It must first keep only the top-K candidates by rigorous Hilbert residual gain.

Only that shortlist may be screened by the known initial reduced electronic density.

This rule prevents the benchmark from turning into an unconstrained fit to one
observable.

## 5. No final-time information in basis construction

The initial basis builder may use:

- the exact initial target wavefunction;
- its initial reduced density;
- the deterministic candidate dictionary.

It may not use:

- exact final populations;
- exact final density matrix;
- final-time trajectory information.

The release benchmark computes the exact final target only after the initial basis has
been selected.

## 6. TDSE defect

The instantaneous dynamical residual is

$$
\mathcal R
=
i\dot\Psi_G-\hat H\Psi_G.
$$

The suite checks:

- finite nonzero defect norm for an incomplete basis;
- finite normalized defect;
- the projection of the defect back onto the current represented Gaussian space is
  tiny relative to the total defect.

This verifies that the reported quantity is predominantly an out-of-span Galerkin
defect rather than a failure of the projected coefficient equation.

## 7. Defect candidate capture

For an orthogonalized new Gaussian pair,

$$
\Delta_c^{\mathrm{TDSE}}
=
\frac{
\sum_a
|\langle g_c^\perp|\mathcal R_a\rangle|^2
}{
\|g_c^\perp\|^2
}.
$$

A test requires this quantity to be positive for a nontrivial admissible candidate.

## 8. Zero-coefficient enrichment

The defect-selected Gaussian pair enters with two exactly zero electronic
coefficients.

The suite verifies:

$$
\Psi_{\mathrm{after}}=\Psi_{\mathrm{before}}
$$

to numerical precision at insertion.

The defect must nevertheless decrease after recomputing the Galerkin coefficient
derivative in the enlarged basis.

## 9. Predicted versus actual defect reduction

The release acceptance compares

$$
\Delta_{\mathrm{predicted}}^{\mathrm{TDSE}}
$$

with

$$
\|\mathcal R_{\mathrm{before}}\|^2
-
\|\mathcal R_{\mathrm{after}}\|^2.
$$

Their relative discrepancy must be below the configured threshold.

## 10. Representation-consistent dynamics

The v0.13 reference initial state is propagated both by:

1. exact 2D TDSE;
2. spinor-complete exact-LVC Gaussian dynamics.

Both start from the identical projected initial wavefunction.

The final reduced-density discrepancy is the projected-state dynamics error.

This remains a separate metric from the error relative to the original intended
coordinate-dependent target state.

## 11. Release thresholds

The reference must satisfy:

```text
monotone residual refinement                     PASS required
initial reduced-density error              <= 0.033
projected-state dynamics density error      <= 2e-4
original-target full-density error          <= 0.033
original-target population L2 error         <= 0.03
coherence phase error                       <= 0.003 rad
generalized norm drift                      <= 1e-4
reference overlap condition number          <= 5e3
defect squared reduction                    >= 1e-8
defect-gain prediction relative error       <= 5e-3
```

These thresholds are regression criteria for this analytic benchmark and are not
universal chemical-accuracy standards.

## 12. PySCF scope

The residual-selection algorithms are demonstrated on the analytic LVC benchmark,
where the target wavefunction can be evaluated explicitly on a two-dimensional grid.

The inherited PySCF backend remains regression-tested, but v0.13 does not claim that a
full molecular TDSE residual is directly available from PySCF.

`V13_PYSCF_RESIDUAL_BRIDGE.md` describes the controlled molecular analogue.
