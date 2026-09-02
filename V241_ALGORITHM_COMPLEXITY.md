# v0.24.1 algorithm complexity

Let `N_AO` be the AO dimension, `N_MO` the scalar MO dimension, `N_det` a CI-vector
dimension, `R` the number of spin-free roots, and
`M=sum_I(2S_I+1)` the complete spin-microstate dimension.

- The v0.24.1 in-core `int2e_p1vxp1` tensor stores `O(N_AO^4)` values and its three
  SOMF contractions cost `O(N_AO^4)`. This is the dominant static integral step and
  intentionally limits the first implementation to small validation molecules.
- AO-to-MO transformation of the three effective one-electron tensors costs
  `O(N_AO N_MO (N_AO+N_MO))` with optimized contractions.
- State-average active-space 1-RDM construction is `R` PySCF RDM evaluations.
- Reduced transition densities require up to `R^2` transition-RDM evaluations;
  spin-ladder work is linear in the number of required ladder steps and orbitals per
  affected CI vector.
- Contracting SOC integrals with all reduced densities costs `O(R^2 N_MO^2)`.
- Complete state-interaction assembly costs `O(M^2)` once reduced amplitudes exist;
  diagonalization and time-reversal audits cost `O(M^3)`.

The pinned OH case has `N_AO=N_MO=6`, `R=3`, and `M=6`, so CI/integral setup dominates
the tiny state-interaction diagonalization. No derivative displacement factor appears
because v0.24.1 is fixed-geometry only. A future trajectory-ready implementation must
separately account for `3A` nuclear coordinates, derivative step ladders, and
cross-geometry overlap/state-tracking work.
