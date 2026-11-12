# Recovered development checkpoint

These documents preserve the September 17, 2026 v0.30.0 checkpoint. Their test
counts and environment descriptions are historical observations, not validation of
the current branch. The archive is not a release.

Source archive: `Gaussian-Nonadiabatic-Dynamics-v0.30.0-development-checkpoint.zip`.
SHA-256: `d15d0b48288514cde3ea78be25bf0f47c4afe78eae11f2787f382ee55cddd1d9`.

The v0.30.1 branch restores the eight missing package modules, four test modules,
and two recomputation examples. It retains current GitHub packaging, stable version,
historical documentation paths, and existing top-level dataset modules. The recovery
manifest in `results/v0301_recovery_manifest.json` records the source hashes.

The recovered `gaussian_dynamics/pre_ml_contract_v294.py` and GitHub's top-level
`v294_preml_contract.py` are different schemas. They are kept under their original
module names; no automatic conversion or v0.29.5 campaign certification is claimed.

The recovered geometry code has one integration correction: coincident endpoints
represent a constant path with identity transport. A dogleg no longer perturbs the
endpoint to fabricate a nonzero path. Other recovered source modules are unchanged.

See [the current iteration](../../../V301_MOLECULAR_HAMILTONIAN.md) for the
new implementation and remaining acceptance work.
