# SOC-Contract Freeze in v0.21.3

v0.21.3 is the last spin-neutral procedural release before the planned analytic SOC
milestone. It freezes what “SOC-ready” means without claiming that SOC physics has been
implemented.

## Frozen contracts

| Area | v0.21.3 rule |
|---|---|
| Model space | Ordered, fixed dimension; no state deletion inside a Gaussian block |
| Multiplets | Optional completeness rule requires exactly $2S+1$ unique components |
| Units | Internal energies in hartree and generalized coordinates in bohr |
| Hamiltonian | Complex Hermitian H with explicit structural residual |
| Physical derivative | Complex Hermitian $K_a=\langle\Phi|\partial_a\hat H|\Phi\rangle$ |
| Connection | Complex anti-Hermitian D; exactly zero in a declared fixed frame |
| Frame | Explicit unitary/isometric matrices; no inferred state-index convention |
| Provenance | Canonical method/model/parameter payload with SHA-256 fingerprint |
| Weak block | Transport last valid density or use zero force if genuinely unseeded |
| Corrector trials | Roll back guide state; commit only the accepted endpoint |
| Initialization | Arbitrary $d$ and $s$ with explicit quadrature weights and frame |
| Cache | Complex data, mandatory provenance key, fixed-frame semantics only |
| Spin-free mode | Permanent and required to remain exactly recoverable |

## Procedure for the first physical SOC model

The v0.22 analytic implementation should not be accepted until all of the following are
true:

1. The basis convention and every state/component label are declared before dynamics.
2. Complete multiplets are included wherever the model claims them.
3. $H_{\mathrm{spin\text{-}free}}$, $H_{\mathrm{SOC}}$,
   $K_{a,\rm spin\text{-}free}$, and $K_{a,\rm SOC}$ are separately testable.
4. The zero-SOC branch reproduces v0.21.3 spin-free H, K, forces, and propagation.
5. H, every K, D, and every frame pass explicit structural tolerances.
6. Constant and coordinate-dependent complex-gauge covariance pass, including an exact
   degeneracy and a vanishing local electronic block.
7. Initial states are projected in a declared global frame or transformed by supplied
   local frame matrices; state indices are never treated as universal identities.
8. Cache/restart identity includes the model space, SOC method, scalar-relativistic
   convention, derivative method, numerical parameters, units, and code/model version.
9. An independent exact-grid reference tests populations, coherences, nuclear density,
   and full-wavefunction fidelity where available.
10. Timestep, Gaussian-basis, electronic-model-space, and sparse-threshold convergence
    are reported separately.

## Procedure for a later ab-initio SOC backend

An ab-initio claim has additional empirical requirements:

- record package and method versions, active space, state averaging, charge,
  multiplicity/component convention, basis/ECP, and scalar-relativistic treatment;
- verify SOC matrix-element units and phase/frame conventions on an independent small
  system;
- provide or rigorously approximate physical SOC derivative operators and identify the
  approximation in provenance;
- validate many-electron cross-geometry overlaps and full-model-space tracking;
- test restart/cache reproduction in the actual electronic-structure runtime;
- separate backend accuracy from Gaussian-basis and time-integration error.

PySCF is not installed or runtime-validated in the v0.21.3 build environment. That fact
is recorded as a boundary, not converted into an inferred validation claim.

## Release decision rule

The first physical SOC release must fail closed: missing derivative terms, incomplete
multiplets, ambiguous units, structural residual violations, provenance mismatch,
moving-frame data in the fixed-frame cache, or a gauge-dependent degeneracy result are
release blockers.

Passing v0.21.3 means the interfaces and procedures are ready for that work. It does not
mean a physical SOC calculation has already been performed.
