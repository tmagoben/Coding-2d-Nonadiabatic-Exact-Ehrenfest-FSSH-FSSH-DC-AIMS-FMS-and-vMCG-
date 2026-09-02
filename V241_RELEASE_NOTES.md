# v0.24.1 release notes

Release date: 2026-08-24

v0.24.1 is a corrective static-SOC release. It implements direct molecular
state-interaction spin-orbit elements with the already pinned PySCF 2.13.1 runtime,
while keeping every moving-nuclei admission gate closed until the required derivative
and cross-geometry evidence exists.

## Added

- `PySCFStateInteractionSOCProviderV241`, using converged common-orbital SA-CASSCF
  roots and returning direct complex `H_spin_free`, `H_soc`, and `H_total` matrices.
- All-electron Breit-Pauli one-electron SOC from `int1e_prinvxp`.
- Two-electron spin-orbit mean field from `int2e_p1vxp1` and the state-averaged
  spin-free 1-RDM, with the explicit `J - 3/2 K_left - 3/2 K_right` contraction.
- Independent finite-sum Clebsch--Gordan evaluation using integer twice-quantum
  numbers; SymPy is not a runtime dependency.
- Wigner--Eckart reconstruction of every component in complete
  `|root,S,M_S>` multiplets.
- PySCF determinant spin-ladder handling when a reference `q=0` Clebsch--Gordan
  coefficient vanishes, including mixed singlet/triplet coverage.
- Explicit time-reversal matrices, root projectors, Kramers-pair diagnostics,
  state order, unit, operator, prefactor, basis, active-space, and runtime identity.
- A pinned OH three-doublet SA-CASSCF(5e,4o)/STO-3G runtime case and an independent
  PySCF JK-driver cross-check of the two-electron SOMF contraction.
- 59 new gates for 315 cumulative release gates.

## Corrected during implementation

PySCF transition RDMs and the SOC integral tensors use an aligned `p,q` contraction
convention. The initial development probe passed zero-imaginary complex CI arrays to
real direct-spin helper routines and transposed the contraction, which understated
the OH matrix. v0.24.1 preserves real scalar CI arrays through ladder/RDM evaluation,
freezes the aligned index convention, and validates the resulting direct matrix by
Hermiticity, time reversal, Kramers pairing, and an independent JK contraction path.

## Claim boundary

The release validates an implementation and a real fixed-geometry molecular
calculation. The OH/STO-3G case is not a basis/method convergence study or an
accuracy benchmark. Physical SOC derivatives, analytic spin-free gradients in this
SOC provider, derivative connections, cross-geometry many-electron overlaps, and
trajectory state tracking are not implemented here. Therefore:

- static PySCF BP-SOMF SOC: validated;
- doublet/Kramers support: validated;
- mixed-multiplicity spin algebra: validated;
- trajectory-ready/live molecular-SOC backend admission: false;
- general ab-initio SOC accuracy: not validated;
- OpenMolcas external admission: still false.

Prism is not copied, imported, packaged, or required. The shipped implementation uses
PySCF's Apache-licensed public APIs and project-native spin algebra.

Recompute with:

```bash
python examples/127_recompute_v0241_pyscf_static_soc.py
python examples/128_recompute_v0241_campaign.py
```
