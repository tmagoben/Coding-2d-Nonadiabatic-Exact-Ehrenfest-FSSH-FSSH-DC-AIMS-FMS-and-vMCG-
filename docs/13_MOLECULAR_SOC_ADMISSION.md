# Molecular SOC backend admission

v0.23.0 turns “supports SOC” into a testable statement with two distinct meanings.

## Protocol-valid versus physically admitted

A provider is protocol-valid when it can supply the full moving-nuclear electronic
contract: spin-free and SOC Hamiltonians, their physical derivatives, derivative
connections, cross-geometry overlaps, convergence state, symmetry data, and immutable
provenance. This is enough to test framework plumbing with deterministic fixtures.

A real backend is admitted only when the protocol passes **and** the source is a
traceable live or externally captured ab-initio calculation with independent reference,
basis, method, rigid-frame, and state-tracking evidence. The two outcomes are separate
fields in every admission report.

## Why static SOC is a separate tier

A static SOC matrix can support a single-geometry calculation. It cannot determine
$K_{a,{\mathrm{SOC}}}$, $D_a$, or cross-geometry state identity. The `static_soc` tier is
therefore rejected for moving-nuclear propagation. `trajectory_ready` is derived only
when every required differential and overlap capability is present.

## Why deterministic replay is part of admission

Electronic-structure installations are expensive and environment sensitive. A replay
captures exact coordinates, component operators, mass matrices, all pair overlaps,
per-record convergence, time reversal, projectors, and complete provenance. The
canonical JSON and deterministic NPZ encodings make accidental or deliberate changes
detectable by SHA-256.

Replay is evidence transport, not evidence creation. Unknown coordinates cannot be
interpolated, and a replay marked `validation_fixture` cannot become a real backend by
passing numerical protocol tests.

## PySCF boundary

The framework probes whether PySCF is importable but does not infer SOC capabilities
from importability. A live bridge must be supplied by a method-specific implementation,
must declare the exact imported version, must expose the complete provider contract,
and must report affirmative SCF, correlated, and SOC convergence for every snapshot.

PySCF was unavailable in the v0.23.0 release environment. No live calculation was run
and no PySCF backend was admitted.

## Next admission milestone

v0.23.1 may admit the first real backend only after all of the following are archived:

1. a method-specific implementation and reproducible environment;
2. explicit molecule, isotope, geometry, charge, state-space, and operator definitions;
3. converged component derivatives and cross-geometry overlaps;
4. independent SOC reference comparison;
5. basis and method convergence ladders;
6. translational and rotational invariance residuals;
7. quantitative state-tracking overlap and assignment-margin evidence;
8. deterministic replay plus a successful real-admission report.

Until then, the analytic even- and odd-electron fixtures remain validation oracles and
spin-free dynamics remains a permanent supported mode.
