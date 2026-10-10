# v0.28.0 Development Validation

The moving-frame layer passes **42/42 focused tests** and its deterministic evidence
campaign passes **50/50 scientific gates**.

Coverage includes frame and transporter unitarity, anti-Hermitian connections,
finite-difference connection/curvature checks, fixed-reference roundtrips, physical
wavefunction equivalence, constant-gauge covariance, transformed TDVP velocities,
implicit-midpoint endpoints, controlled lifecycle events, and an independent gauge-link
lattice Hamiltonian with Hermiticity, unitary-similarity, action, and propagation checks.

Negative controls reject non-flat curvature, missing trivialization, nonunitary gauges,
non-Hermitian generators, invalid phase Hessians, and unsupported lattice mass tensors.

The evidence intentionally keeps these claims false: nonzero-curvature connections,
live molecular-SOC trajectories, general ab-initio SOC-dynamics accuracy, full AIMS
branching semantics, and release-level certification of v0.28.0. v0.27.0 remains the
sealed release until the complete inherited campaign is rerun successfully.
