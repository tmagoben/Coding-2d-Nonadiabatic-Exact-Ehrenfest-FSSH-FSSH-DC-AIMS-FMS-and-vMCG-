# v0.30.0 Molecular Electronic-Frame Bridge

## Classification

v0.30.0 is a **development checkpoint**. It begins the molecular integration
promised after v0.29.4, but it does not replace the certified v0.27.0 release and
does not yet claim a complete live ab-initio trajectory.

## Purpose

The checkpoint binds the v0.29.3 coordinate-dependent electronic geometry to an
existing molecular electronic-structure provider. At every generalized nuclear
coordinate `q`, the provider must return

\[
H(q),\qquad K_\mu(q)=\langle\Phi|\partial_\mu\hat H|\Phi\rangle,
\qquad D_\mu(q)=\langle\Phi|\partial_\mu\Phi\rangle,
\qquad M(q),
\]

with `H` and every physical operator derivative `K_mu` Hermitian,
`D_mu` anti-Hermitian, and `M` real symmetric positive definite. `H` is the total
electronic Hamiltonian, including SOC whenever the provider provenance enables it.

The bridge preserves the distinction between the physical derivative operator
`K_mu` and the derivative of a matrix written in a moving frame. The latter contains
connection commutators and must never be substituted for `K_mu`.

## Two independent transport paths

For a path from `a` to `b`, the NAC/connection route integrates

\[
\frac{dW}{ds}=-\dot R^\mu D_\mu(R(s))W,
\qquad W(0)=I.
\]

The overlap route evaluates the physical finite-manifold matrix

\[
O_{ba}=\langle\Phi(b)|\Phi(a)\rangle
\]

and extracts its unitary polar factor only after contraction, retention,
conditioning, principal-angle, unitarity, positivity, and reciprocity checks. The raw
overlap is retained as evidence but is never used directly to transport coefficients.

For sufficiently local segments, the two unitary transports must agree within the
configured audit tolerance. A disagreement fails closed because it signals a NAC
orientation, phase, state-order, truncation, or provider-consistency defect.

## Molecular context for the spawning dataset

The v0.29.4 exact residual oracle remains unchanged and authoritative. v0.30.0 may
attach a separate molecular context record containing

- provider and bridge fingerprints;
- stable packet IDs and point fingerprints;
- sorted adiabatic energy spectra;
- minimum adiabatic gaps;
- total physical operator-gradient norms;
- SOC provenance and the inherited trajectory-group split.

Raw Hamiltonians, derivative connections, SOC components, electronic coefficients,
and gauge-dependent vectors are not model features. This preserves the initial safety
boundary

\[
\boxed{\text{ML shortlist/rank}\;\longrightarrow\;\text{exact spawning oracle}}.
\]

## Fail-closed boundaries

The bridge rejects providers that change coordinate or electronic dimension, return a
point at the wrong coordinate, violate Hermiticity/anti-Hermiticity/mass positivity,
omit provenance, supply spectrally expansive physical overlaps, violate reciprocal
overlap/transport identities, or fail a required v0.23.3 numerical-identity gate. A
production numerical identity must be declared by the actual outer provider; it is not
silently inherited through a gauge, tracking, or interpolation wrapper.

## Explicit non-claims

This checkpoint does not yet provide

- a full nuclear-kinetic plus electronic molecular `H psi` trajectory loop;
- patch-aware implicit-midpoint propagation across changing atlas membership;
- an automatic conical-intersection atlas or path optimizer;
- production live PySCF NAC/SOC trajectories in the curved basis;
- molecular accuracy claims;
- ML-controlled spawning or unrestricted AIMS branching.
