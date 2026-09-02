# v0.24.1 validation

## Validated implementation properties

- Exact PySCF 2.13.1 distribution and module identity.
- Availability of `int1e_prinvxp`, `int2e_p1vxp1`, spin-separated transition 1-RDM,
  and determinant creation/annihilation APIs.
- Converged OH ROHF and equal-weight three-root SA-CASSCF(5e,4o)/STO-3G.
- Direct nonzero six-state complex `H_soc`, not an eigenpair reconstruction.
- Explicit rank-five SOMF contraction versus a separate PySCF JK-driver path.
- One-electron/two-electron AO antisymmetry and effective-MO Hermiticity.
- Exact `H_total=H_spin_free+H_soc`, unitary eigensystem, and reconstruction.
- Complete doublet multiplets, odd-electron time-reversal square, and Kramers pairs.
- Mixed singlet/triplet Wigner reduction, including the zero-`q=0`-CG spin-ladder
  path, plus doublet/quartet completeness and rank-one spin selection rules.
- Native-Boolean static-only capability and negative admission controls.

The canonical OH evidence reports a SOC Frobenius norm near 139 cm^-1 and a largest
absolute matrix element near 51.4 cm^-1. These values identify the smoke calculation;
they are not experimental or complete-basis reference claims.

## Deliberately unvalidated

- analytic or finite-difference physical `dH_soc/dR` for production;
- SOC-aware analytic spin-free gradients in this provider;
- derivative connections and the frozen NAC convention at moving geometries;
- cross-geometry many-electron overlaps and connected state tracking;
- basis-set, active-space, dynamic-correlation, scalar-relativistic, or operator-model
  convergence;
- agreement with an independent production quantum-chemistry implementation;
- general spectroscopic or ab-initio SOC accuracy;
- trajectory-ready/live molecular-SOC backend admission.

The v0.24.1 release campaign contains 315 gates: 256 inherited v0.24.0 gates, 39 real
PySCF static-SOC runtime gates, and 20 corrective spin-algebra/admission gates.
