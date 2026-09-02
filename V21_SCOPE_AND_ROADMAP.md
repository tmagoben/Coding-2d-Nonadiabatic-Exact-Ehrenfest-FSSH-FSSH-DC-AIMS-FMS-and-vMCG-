# Scope and roadmap

## v0.21 capabilities

v0.21 supports:

- complex Hermitian electronic Hamiltonians;
- complex Hermitian physical Hamiltonian-derivative operators;
- complex anti-Hermitian derivative connections;
- arbitrary electronic block dimension;
- smooth local \(U(s)\) gauge transformations;
- full-subspace Procrustes alignment;
- Wilson-loop diagnostics;
- block-sparse \(S/H/T\);
- dynamic sparse topology;
- sparse moving-basis propagation;
- molecular provider adapters.

## v0.21 exclusions

v0.21 does not yet implement:

- an SOC Hamiltonian;
- spin labels as a special propagation concept;
- production AIMS matrix elements;
- a real PySCF v0.21 runtime trajectory;
- production asynchronous electronic-structure scheduling.

## Why an SOC-neutral core is preferable

A future calculation may use

$$
H_e
=
H_{\rm spin-free}
+
H_{\rm SOC}.
$$

A spin-free calculation simply has

$$
H_{\rm SOC}=0.
$$

The Gaussian engine should not care which physical terms produced the final Hermitian
electronic operator.

The intended structure is therefore:

```text
same Gaussian engine
same block sparse graph
same gauge algebra
same convergence machinery
different optional electronic operator backends
```

## Recommended next milestone

Before introducing physical spin terms, the strongest remaining validation target is a
small real molecular calculation using the v0.21 interfaces:

1. SA-CASSCF energies;
2. gradients;
3. NACs;
4. many-electron cross-geometry overlaps;
5. conversion to full complex operator matrices;
6. sparse threshold convergence;
7. timestep convergence;
8. cache reproducibility.

## Future SOC route

When SOC is introduced, it should enter at the electronic layer:

$$
H(q)
=
H_0(q)
+
H_{\rm SOC}(q).
$$

The first SOC benchmark should be analytic and exactly reproducible.

Only after that passes should an ab-initio SOC backend be connected.

## Permanent non-SOC mode

Ordinary spin-free nonadiabatic dynamics should remain a first-class mode indefinitely.

SOC is an optional physical model, not the definition of the framework.
