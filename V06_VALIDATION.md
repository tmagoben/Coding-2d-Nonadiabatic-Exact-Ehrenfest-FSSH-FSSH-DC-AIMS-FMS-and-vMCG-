# v0.6 Validation Contract

## State assignment

- known root permutation is recovered;
- known sign inversion is removed;
- transformed NAC tensor follows
  $d'_{ij}=p_i^*p_jd_{\pi(i)\pi(j)}$;
- ambiguous equal-score assignments are flagged.

## Many-electron CASSCF overlap

- active CI vectors are embedded with doubly occupied core orbitals;
- nonorthogonal determinant overlap is evaluated in the full core+active orbital
  space;
- a test with explicit core-active orbital rotation reproduces the exact
  many-electron overlap;
- a synthetic root swap/sign flip is recovered directly from the many-electron
  overlap matrix.

## Subspace/gauge transport

- unitary Procrustes alignment exactly removes a known unitary basis rotation;
- principal-overlap singular values identify an unchanged selected subspace;
- overlap unitarity defect is reported.

## Overlap-derived NAC

- a known two-state rotation produces the expected anti-Hermitian generator as the
  geometry step tends to zero.

## Tracked PySCF backend

- first geometry initializes labels without arbitrary reassignment;
- second geometry can reorder energies/gradients/NACs by tracked identity;
- ambiguous tracking raises by default;
- `reset_tracking()` clears electronic-state history.

## Tracked scan

- a sequentially produced tracked 1D scan can be queried afterward in arbitrary order;
- interpolation preserves exact NAC antisymmetry.

## Full regression

All earlier v0.1-v0.5 tests must continue to pass.

## Real PySCF requirement

Before research use, a machine with PySCF installed must additionally check:

1. overlap continuity versus geometry-step refinement;
2. overlap-derived directional NAC versus analytic PySCF NAC projection;
3. state identities against orbital/configuration character;
4. tracking around any closed loop relevant to geometric phase;
5. stability to active-space and state-manifold enlargement.
