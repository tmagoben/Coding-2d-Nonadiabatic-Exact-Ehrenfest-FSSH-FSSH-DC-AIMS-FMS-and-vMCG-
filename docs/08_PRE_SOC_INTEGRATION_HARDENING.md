# Pre-SOC integration hardening in v0.21.2

v0.21 proved the complex block/gauge algebra. v0.21.2 connects the pieces that must work
before a physical SOC term is introduced.

The key design decision is that **SOC will be optional electronic physics**, not a new
propagation architecture.

The release therefore hardens six spin-neutral interfaces:

1. unequal-width Gaussian pair algebra;
2. coefficient-coupled nuclear guidance;
3. complete electronic-block birth/prune operations;
4. physical electronic observable operators;
5. full-subspace continuity diagnostics;
6. complex dtype preservation.

See `V212_DERIVATIONS.md` for the step-by-step mathematics and
`V212_VALIDATION.md` for numerical acceptance values.

## Boundary found by the v0.21.3 audit

This document records the v0.21.2 state. A subsequent exact-degeneracy audit found that
the low-amplitude lowest-eigenvector fallback was not representation covariant, and a
strict-invariant audit found that the old `numpy.allclose(..., atol=...)` calls inherited
a relative tolerance that could admit a materially non-Hermitian operator.

v0.21.3 closes both defects and freezes the surrounding operator/model-space,
initialization, and cache procedures. See `09_SOC_CONTRACT_FREEZE.md`.

A later physical SOC release may then add

$$
H(q)=H_0(q)+H_{SOC}(q)
$$

through the frozen complex electronic-operator interface, but only with both the
Hamiltonian and physical derivative terms plus explicit provenance.

A real ab-initio SOC claim should still wait for a separate real electronic-structure
runtime validation, but an exactly reproducible analytic SOC benchmark no longer needs
to wait on that empirical backend milestone.
