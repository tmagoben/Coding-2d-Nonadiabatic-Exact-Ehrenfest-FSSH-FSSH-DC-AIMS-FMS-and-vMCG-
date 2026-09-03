# v0.21.4 Release Notes

v0.21.4 is the **differential-provider and deterministic-restart certification
release**. It remains spin-physics neutral and preserves every frozen v0.21.3
electronic contract.

This release deliberately contains **no physical SOC Hamiltonian or SOC derivative**.

## Added

- a centered cross-geometry provider audit for physical Hamiltonian derivatives K,
  derivative connections D, overlap isometry, pointwise structure, and provenance;
- adversarial fixtures whose matrices are pointwise valid but differentially wrong;
- versioned self-consistent block checkpoints with strict provider and settings
  fingerprints;
- a SHA-256 digest over the canonical manifest and every numerical array;
- checkpointed Gaussian UIDs, widths, positions, momenta, coefficients, density-guide
  masks/densities/counters, and sparse active UID edges;
- deterministic dense, sparse, moving-complex-frame, and adaptive-lifecycle restart;
- an explicit zero-SOC wrapper that routes exact complex zero matrices through the
  frozen H/K composition contract;
- a canonical 21-gate v0.21.4 campaign inheriting all v0.21.3 acceptance checks.

## Restart semantics

A restart is accepted only when checkpoint format, integrity digest, provider
provenance fingerprint, propagation-settings fingerprint, model dimension, coordinate
dimension, dt, and array structure all agree. Provider wavefunction objects are not
serialized; the declared provider reconstructs snapshots at checkpoint geometries.

Adaptive policies receive the global step after restart. Sparse hysteresis is restored
by stable Gaussian UID edges, not transient array indices. Density guidance resumes
with its accepted densities and diagnostic counters, including an exactly zero local
coefficient block with a valid retained guide density.

## Compatibility

The v0.21.3 model-space, provenance, H/K/D, density-guidance, initialization, caching,
spin-free propagation, and acceptance contracts remain intact. Existing v0.21.2 and
v0.21.3 entry points are preserved.

## Deferred to v0.22 and later

- a physical analytic H_SOC and K_SOC;
- an independently propagated exact-grid SOC reference;
- ab-initio SOC matrices or derivatives;
- a real PySCF runtime claim;
- production AIMS equations or asynchronous scheduling.

Use this release description:

> Pre-SOC, complex representation-neutral block-sparse Gaussian nonadiabatic dynamics
> with cross-geometry provider certification, deterministic provenance-checked restart,
> degeneracy-safe density guidance, and exact zero-SOC integration rehearsal.

Do not describe v0.21.4 as SOC dynamics or production AIMS.

