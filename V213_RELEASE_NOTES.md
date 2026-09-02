# v0.21.3 Release Notes

v0.21.3 is the **SOC-contract-freeze release**. It removes the last known
representation-dependent fallback and makes the electronic interface auditable before
any physical spin-orbit term is added.

This release deliberately contains **no physical SOC Hamiltonian**.

## Added

- explicit scaled Frobenius residuals for Hermiticity, anti-Hermiticity, isometry, and
  mass symmetry;
- fixed-dimension electronic model-space descriptors with optional complete-multiplet
  enforcement;
- canonical electronic-operator provenance and SHA-256 fingerprints;
- explicit hartree/bohr internal-unit rules and cm$^{-1}$ conversion helpers;
- a composition contract for
  $H=H_{\rm spin\text{-}free}+H_{\rm SOC}$ and
  $K_a=K_{a,\rm spin\text{-}free}+K_{a,\rm SOC}$;
- transported density-matrix guidance for exact degeneracies and weak local blocks;
- parent-density inheritance for a newly inserted zero electronic block;
- arbitrary-state, arbitrary-nuclear-dimension fixed-frame grid projection;
- a fixed-frame complex operator cache keyed by the complete provenance fingerprint;
- transactional density-guide checkpoints across predictor/corrector trials;
- a canonical 20-gate v0.21.3 acceptance campaign.

## Changed

- `ElectronicOperatorPointV21.validate` and snapshot/frame validation no longer depend
  on NumPy's implicit `allclose` relative tolerance.
- `dH_dq` is retained for compatibility, while
  `hamiltonian_derivative_operator_q` states its physical meaning explicitly.
- v0.21.3 self-consistent propagation uses density-matrix guidance by default.
- fixed-frame cache entries contain complex H/K/D data, mass matrices, state vectors,
  metadata, and a mandatory provider fingerprint.

## Retired

The low-amplitude `lowest_eigenvector` guidance policy is rejected. At an exact
degeneracy it selected an arbitrary basis vector and produced gauge-dependent nuclear
forces. A weak block now uses its transported last valid density, or zero force if it
has never had a physical or inherited guide density.

## Compatibility

The Gaussian block propagation architecture, sparse graph, spin-free provider path,
and v0.21.2 release campaign remain intact. Spin-free dynamics is a permanent mode, not
a temporary special case.

## Deferred to v0.22 and later

- a physical analytic SOC Hamiltonian and its physical derivative operators;
- an exact-grid SOC dynamics reference;
- an ab-initio SOC backend;
- a real PySCF SOC runtime claim;
- production AIMS nuclear equations or production asynchronous scheduling.

Use this release description:

> SOC-contract-frozen, complex representation-neutral block-sparse Gaussian
> nonadiabatic dynamics with strict electronic invariants, degeneracy-safe density
> guidance, fixed model spaces, arbitrary-state initialization, and provenance-safe
> caching.

Do not describe v0.21.3 as SOC dynamics or production AIMS.
