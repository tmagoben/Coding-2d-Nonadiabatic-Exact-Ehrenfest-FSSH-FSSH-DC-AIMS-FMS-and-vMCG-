# Restricted time-dependent-variational SOC dynamics

v0.25.0 adds a symmetric trajectory layer for one constant-mass canonical nuclear
packet coupled to a complete SOC spinor. The implementation is a controlled TDVP
restriction: it is not the full coupled multi-Gaussian variational method.

## Algorithm

For each signed step, the runner evaluates the start force
`F_a=-c^dagger K_a c`, applies a nuclear half kick, drifts with `M^-1`, evaluates the
endpoint operator snapshot, and constructs the cross-geometry overlap. If
`O=U Sigma V^dagger`, the endpoint amplitudes are

$$
c_1=e^{-iH_1h/2}V U^\dagger e^{-iH_0h/2}c_0.
$$

Since `V U^dagger = W^dagger`, this is endpoint Strang propagation with the unitary
right-to-left polar factor reversed into start-to-end transport. The endpoint
electronic force completes the second nuclear half kick.

## Why SVD-polar

The raw overlap is a physical contraction and is retained. Its singular values test
manifold retention and conditioning. The polar factor is the unitary object that
transports amplitudes. SVD computes both in one robust factorization, so the project
uses polar transport *via* SVD rather than choosing one over the other.

## Scope

Velocity Verlet is admitted only for the released constant-mass canonical variables.
The future full multi-Gaussian TDVP must use an implicit midpoint/discrete variational
solve with metric, gauge-null-space, constraint, and nonlinear-convergence controls.
Real PySCF molecular SOC trajectories remain closed until an arbitrary-geometry
provider supplies all Cartesian `H`, physical `K`, `D`, mass, overlaps, and independent
accuracy evidence.

See `V250_VARIATIONAL_SOC_DYNAMICS.md`, `V250_INTEGRATOR_DECISION.md`, and
`V250_PROGRAM_ARCHITECTURE.md` for the complete equations and data flow.
