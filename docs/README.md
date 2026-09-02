# Documentation map

The documentation is organized by concept rather than only by release number.

1. `01_MATHEMATICAL_FOUNDATIONS.md` — step-by-step derivation from the molecular TDSE to the moving nonorthogonal Gaussian equations.
2. `02_COMPLEX_GAUGE_FRAMEWORK.md` — complex U(s) representation covariance, physical operator derivatives, subspace tracking, and Wilson loops.
3. `03_BLOCK_SPARSE_ALGORITHM.md` — derivation of block S/H/T, the gauge-invariant sparse score, and moving-basis propagation.
4. `04_ALGORITHMIC_COMPLEXITY.md` — time and memory complexity in N, s, E, M, d, and electronic-structure cost.
5. `05_VALIDATION_STRATEGY.md` — layered acceptance criteria and failure isolation.
6. `06_MOLECULAR_BACKENDS.md` — molecular electronic-structure interfaces, PySCF conventions, caching, and tracking boundaries.
7. `07_SCOPE_AND_ROADMAP.md` — validated real spin-free runtime scope and the next method-specific SOC milestone.
8. `08_PRE_SOC_INTEGRATION_HARDENING.md` — v0.21.2 integration layers and the defect later closed by v0.21.3.
9. `09_SOC_CONTRACT_FREEZE.md` — the v0.21.3 operator, model-space, guidance, initialization, cache, and release gates required before physical SOC.
10. `10_DIFFERENTIAL_AND_RESTART_CERTIFICATION.md` — the v0.21.4 cross-geometry provider, deterministic restart, and zero-SOC rehearsal procedures.
11. `11_PHYSICAL_ANALYTIC_SOC.md` — the v0.22.0 singlet–triplet and Kramers-doublet physical-SOC contract and claim boundary.
12. `12_CORRECTIVE_SOC_HARDENING.md` — the v0.22.1 derivative, symmetry, exact-grid, and convergence admission barrier.
13. `13_MOLECULAR_SOC_ADMISSION.md` — the v0.23.0 capability, evidence, replay, and real-admission boundary.
14. `14_RAW_EVIDENCE_ADMISSION.md` — the v0.23.1 receipt, raw-observation, manifold-tracking, and executable-validator boundary.
15. `15_PYSCF_RUNTIME_AND_OVERLAP_CORRECTIONS.md` — the v0.23.2 real runtime, NAC orientation, finite-manifold overlap, and admission corrections.
16. `16_TRANSPORT_AND_COMPATIBILITY_HARDENING.md` — the v0.23.3 raw-overlap/unitary-transport split, replay migration, convention identities, complete manifolds, and runtime profiles.
17. `RELEASE_HISTORY.md` — development arc through v0.23.3.
