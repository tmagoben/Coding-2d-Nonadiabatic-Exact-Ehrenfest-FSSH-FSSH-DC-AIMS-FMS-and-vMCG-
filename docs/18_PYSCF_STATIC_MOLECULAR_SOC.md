# PySCF static molecular SOC

v0.24.1 can compute direct fixed-geometry molecular SOC elements from converged
PySCF common-orbital SA-CASSCF roots.

## Reproduce the admitted evidence

Use the hash-locked CPython 3.12/Linux wheel environment described by
`requirements-pyscf-v241-linux-x86_64-py312.txt`, fix BLAS/OpenMP threads to one, and
run:

```bash
python examples/127_recompute_v0241_pyscf_static_soc.py
```

The evidence JSON contains `H_spin_free`, direct `H_soc`, `H_total`, SOC eigenvectors,
eigenvalues, the complete state order, runtime fingerprints, molecular/method identity,
integral summaries, all 39 runtime gates, and explicit false trajectory/accuracy
claims.

## Programmatic use

Given a converged common-orbital PySCF CASSCF object `mc` and the already verified
runtime fingerprint:

```python
provider = PySCFStateInteractionSOCProviderV241(
    mc,
    environment_sha256=runtime_fingerprint.environment_sha256,
    root_labels=("D1", "D2", "D3"),
    root_spin_twice=(1, 1, 1),
    basis_label="STO-3G",
)
result = provider.evaluate_static_soc()
H_soc_hartree = result.matrices.H_soc
state_order = result.matrices.state_order
```

Never discard `state_order`, `convention`, or `identity`: matrix elements are
basis-, phase-, ordering-, operator-, and method-specific. Compare invariant spectra
or explicitly aligned matrices, not untracked raw elements from unrelated runs.

## Why dynamics still refuses this provider

A force or moving-basis propagation needs more than `H_soc(R0)`. It needs the
physical fixed-frame operator derivative, spin-free derivative, electronic
connection, cross-geometry correlated overlaps, and certified state transport at
every required geometry. Returning zeros for those quantities would silently change
the physics. v0.24.1 therefore raises at every moving-geometry method and advertises
only the `static_soc` tier.

The next scientifically meaningful step is transported finite-difference BP-SOMF
evidence over a small connected geometry scan, followed by independent operator and
basis/method convergence checks. Only then should trajectory admission be reconsidered.
