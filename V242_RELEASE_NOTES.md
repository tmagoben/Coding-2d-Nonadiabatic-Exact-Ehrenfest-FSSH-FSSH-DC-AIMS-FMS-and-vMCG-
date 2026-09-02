# v0.24.2 release notes

Release date: 2026-08-24

v0.24.2 is the connected-geometry SOC differential-preview release. It replaces the
v0.24.1 production rank-five two-electron SOMF contraction with PySCF's direct JK
driver, retains overlap-capable SA-CASSCF snapshots at neighboring geometries, and
forms separately transported spin-free and SOC centered differences in a complete
doublet microstate space.

## Added

- A direct-JK implementation of
  `J[D] - 3 K_left[D]/2 - 3 K_right[D]/2` for `int2e_p1vxp1`. The production
  path stores only the three AO output matrices; the explicit rank-five tensor is
  used once at the small OH center as a validation oracle.
- `PySCFSOCGeometrySnapshotV242`, binding geometry, roots, common scalar orbitals,
  CI vectors, state-average density, BP-SOMF integrals, state-interaction matrices,
  spin order, calculation input, and exact runtime identity.
- Exact restricted-CASSCF many-electron cross-geometry overlaps lifted from root
  space into complete `|root,S,M_S>` spaces without mixing different `S` or `M_S`.
- Degenerate-safe right-to-left unitary polar transport. This is essential for the
  first two OH doublet roots, which undergo large arbitrary rotations between
  independently converged SA-CASSCF calculations.
- Separate transported centered differences for `H_spin_free`, `H_soc`, and their
  total. Each reported component is recomputed from and bound to its serialized
  transported endpoint matrices.
- A three-level `0.08, 0.04, 0.02` bohr OH bond-coordinate ladder on the observed
  second-order truncation plateau, with Richardson estimates retained as evidence.
- Six serialized endpoint receipts and fingerprint, environment, state-order,
  state-average, and signed-geometry binding for all three displacement pairs.
- 60 real-runtime gates and 25 adversarial/core gates, giving 85 new and 400
  cumulative release gates.

## Corrected during implementation

Root-by-root sign alignment was rejected for the OH scan after the nominally
degenerate first two roots were observed to rotate strongly even for small geometry
changes. Operator transport now uses the full retained-manifold polar factor. The
anti-Hermitian part of the polar-aligned overlap slope is recorded only as a local
parallel-transport-gauge preview; it is not relabeled as a physical derivative
connection.

The initial differential receipt stored endpoint fingerprints but serialized only
the center snapshot. The final schema retains all six compact endpoint receipts and
verifies that every fingerprint and signed displacement matches its derivative
record. It also binds `K_spin_free` and `K_soc` independently to the corresponding
transported endpoint matrices, preventing equal-and-opposite component tampering.

## Claim boundary

Validated in v0.24.2:

- real PySCF 2.13.1 direct-JK BP-SOMF execution;
- agreement with the explicit small-system SOMF oracle;
- connected OH geometry snapshots and many-electron overlaps;
- complete three-doublet/six-microstate polar transport;
- separate transported OH bond-coordinate spin-free and SOC derivative previews;
- centered-difference second-order behavior over the frozen ladder.

Not validated or admitted:

- a continuous physical derivative connection/NAC from the SOC spinor states;
- the full `3N` Cartesian derivative tensor;
- analytic spin-free or SOC derivatives in this provider;
- a real mixed-multiplicity molecular runtime case;
- basis/method convergence or general ab-initio SOC accuracy;
- a live or trajectory-ready molecular-SOC backend.

Recompute with:

```bash
python examples/129_recompute_v0242_pyscf_differential_soc.py
python examples/130_recompute_v0242_campaign.py
```
