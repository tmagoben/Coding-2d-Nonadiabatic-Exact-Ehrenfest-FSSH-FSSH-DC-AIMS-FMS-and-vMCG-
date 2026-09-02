# Differential Provider and Restart Certification

v0.21.4 adds the final analytic-provider and state-continuity checks before the planned
v0.22 physical SOC milestone. It preserves the v0.21.3 frozen operator contract and
introduces no physical spin term.

## Provider admission procedure

A provider intended for a later physical extension must pass, at representative
geometries:

1. the v0.21.3 model-space, provenance, unit, and structural checks;
2. centered overlap-transported K checks in every nuclear coordinate;
3. centered raw-overlap D checks in every nuclear coordinate;
4. cross-geometry overlap shape, finiteness, and isometry checks;
5. invariant provider provenance across all displaced evaluations;
6. negative controls demonstrating that wrong K and wrong D are detected;
7. a documented displacement-refinement study when the provider is numerical.

Passing at one geometry is not universal validation. Molecular backends should select a
set spanning ordinary regions, avoided crossings, near-degeneracies, and any region
where tracking changes are plausible.

## Restart admission procedure

A deterministic restart must preserve both visible wavefunction data and hidden
controller state. v0.21.4 therefore requires:

- stable Gaussian UIDs and the full q/p/A/C block state;
- accepted density-guide masks, matrices, and counters;
- sparse active UID edges used by hysteresis;
- global step, time, and dt;
- exact checkpoint-format identity;
- provider-provenance and propagation-settings fingerprints;
- a canonical full-state integrity digest;
- reconstruction of provider snapshots rather than serialization of backend-specific
  wavefunction objects;
- uninterrupted/segmented equivalence in fixed, sparse, and moving complex frames;
- insertion/pruning equivalence across a segment boundary.

## Zero-SOC admission procedure

Before a physical optional term is added, its disabled branch must route through the
same composition seam and reproduce the permanent spin-free path. v0.21.4 requires
exact operator equality and near-roundoff trajectory equality with explicit complex
zero H_SOC and K_SOC.

## What remains for v0.22

The first physical SOC release still needs a reproducible analytic H_SOC, its physical
K_SOC derivatives, a complete fixed model space and spin convention, complex-gauge
covariance, an independent exact-grid dynamics reference, and convergence in timestep,
Gaussian basis, and sparse threshold. Ab-initio SOC remains a later validation layer.

