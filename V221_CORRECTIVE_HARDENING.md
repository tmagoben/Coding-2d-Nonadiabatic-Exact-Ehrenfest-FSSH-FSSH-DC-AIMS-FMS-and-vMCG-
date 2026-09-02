# v0.22.1 corrective SOC hardening

## Why this release exists

v0.22.0 established a physical analytic-SOC path, but several procedures were too
specialized or too weak for admitting a molecular backend. v0.22.1 turns each gap into
an explicit fail-closed contract and an adversarial regression.

## Component-resolved derivative certification

For every coordinate `a` and displacement `h`, neighboring operators are transported
to the center electronic frame before forming centered differences:

\[
\delta_a H_X(h)=
\frac{\widetilde H_X(q+h e_a)-\widetilde H_X(q-h e_a)}{2h},
\qquad X\in\{\mathrm{sf},\mathrm{SOC}\}.
\]

The audit separately requires

\[
\delta_aH_{\mathrm{sf}}\approx K_{a,\mathrm{sf}},
\qquad
\delta_aH_{\mathrm{SOC}}\approx K_{a,\mathrm{SOC}}.
\]

This is stronger than checking only their sum or one density-contracted force. The
negative fixture adds `Delta` to one component and subtracts it from the other while
choosing `Delta` orthogonal to the sampled density. The old total-K and scalar-force
checks pass; the new component checks reject it with residuals near `6.69e-4`.

The production audit is dimension neutral: it uses the emitted coordinate length and
state dimension, iterates over every K matrix, and requires no provider `.config`.

## Symmetry admission contract

An SOC provider must now emit a `SOCSymmetryContractV221` containing:

- declared even or odd electron parity;
- the numerical antiunitary time-reversal matrix `J`;
- the numerical physical projector family;
- the external-field declaration.

Admission requires complete multiplets, a single charge and electron-parity sector,
fermionic-declaration agreement, zero external magnetic field,

\[
J^\dagger J=I,\qquad JJ^*=\begin{cases}+I&\text{even},\\-I&\text{odd},\end{cases}
\]

and Hermitian, idempotent, mutually orthogonal projectors that resolve the identity.
The numerical `J` and projectors are part of provenance identity, not only descriptive
labels. Mixed singlet/doublet sectors, a nonunitary `J` that nevertheless has the right
square, and changed projector provenance are independent rejection controls.

## Exact-grid oracle contract

The exact-grid solver is deliberately one-dimensional and fixed-frame. It now proves
those preconditions from provider output. Mass comes from `mass_matrix_q_au`; all grid
points must agree on one positive scalar mass. A moving electronic frame is rejected
because a split operator omitting the associated gauge-connection terms would not be a
valid independent reference.

For a static potential, the potential half-step matrices and kinetic phase are
precomputed once. The public single-step function remains available and is tested
against the precomputed trajectory. Saved times always include `steps * dt`.

## Convergence conditions added

The v0.22.0 Gaussian/grid comparison was retained, and two orthogonal convergence
tests were added:

1. The physical triplet population is recomputed in prescribed 1-, 3-, and 5-Gaussian
   bases. The 3-to-5 change must be below `1e-8`, and its ratio to the 1-to-3 change must
   be below `0.05`.
2. A doublet SOC trajectory is recomputed for four decreasing sparse thresholds. Active
   edge counts must not decrease, coefficient errors against the dense solution must
   not increase beyond `1e-14`, and the finest error must be below `1e-12`.

These are analytic framework gates, not evidence of molecular basis-set convergence.

## Gate to v0.23

A future molecular SOC backend should be admitted only if it satisfies this same
component-resolved differential audit, symmetry/provenance admission, fixed unit and
model-space identity, molecular method/basis convergence, independent reference-data
comparisons, and deterministic restart identity. v0.22.1 itself makes no such backend
claim.

