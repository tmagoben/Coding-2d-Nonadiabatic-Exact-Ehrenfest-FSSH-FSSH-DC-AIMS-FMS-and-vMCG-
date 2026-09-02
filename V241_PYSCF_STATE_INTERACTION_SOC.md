# v0.24.1 PySCF state-interaction SOC

## Operator identity

v0.24.1 implements an all-electron Breit-Pauli spin-orbit mean-field operator in a
real scalar-orbital basis. PySCF supplies Cartesian AO tensors in `(x,y,z)` order:

```text
one-electron: int1e_prinvxp
two-electron: int2e_p1vxp1
```

For a state-averaged spin-free AO density `P`, the two-electron mean-field tensor is

```text
F_2e = J(P) - 3/2 K_left(P) - 3/2 K_right(P).
```

The effective spatial SOC tensor is

```text
h_SOC = -i C^T [ (1 / (2 c^2)) (W_1e - F_2e) ] C,
```

where `C` is the common real SA-CASSCF MO coefficient matrix. The factor
`1/(2c^2)` is included exactly once. No scalar-relativistic correction or ECP-SOC
term is mixed into this operator identity.

## State-interaction basis

Each spin-free root `I` declares `(E_I,S_I,M_ref)`. The output basis contains every
component

```text
|I,S_I,M_S>,  M_S = S_I, S_I-1, ..., -S_I,
```

ordered first by root and then by decreasing `M_S`. This makes each doublet an
adjacent `(+1/2,-1/2)` Kramers block. The complete order is stored in both the result
and the frozen v0.23.3 convention object; consumers must not infer or reorder it.

PySCF's spin-separated transition 1-RDMs provide a nonzero `q=0` reference component.
If its Clebsch--Gordan coefficient vanishes, normalized `S+` or `S-` determinant
operations move both roots to another common `M_S`. The reduced rank-one density is
then expanded with the Wigner--Eckart theorem. This handles integer-spin `M_S=0`
roots without division by zero and supports same-parity mixed multiplicities.

All angular momenta are represented internally as integers `2S` and `2M_S`. A finite
Racah sum evaluates Clebsch--Gordan coefficients, eliminating floating half-integer
ambiguity and any SymPy dependency.

## Direct output contract

The provider returns, rather than reconstructs, these matrices:

- `H_spin_free` in hartree;
- `H_soc` in hartree;
- `H_total = H_spin_free + H_soc`;
- the eigenvalues/eigenvectors of `H_total`;
- the full state order, time-reversal matrix, and root projectors;
- the AO/MO integral identity and state-average density;
- runtime, molecular, method, basis, active-space, unit, and convergence provenance.

The matrix is rejected rather than symmetrized if its raw Hermiticity residual exceeds
the tolerance. Odd-electron results additionally require the antiunitary time-reversal
square `Theta^2=-1` and numerical Kramers pairing.

## Static-only boundary

`evaluate_static_soc()` is the only evaluation entry point admitted in v0.24.1.
`components(q)`, `evaluate_snapshot(q)`, and `snapshot_overlap(left,right)` raise
because the provider has no physical SOC derivatives or cross-geometry wavefunction
overlaps. Its `MolecularSOCCapabilitiesV230` receipt is exactly `static_soc`; the
derived `trajectory_ready` and real-backend-admission properties remain false.

The frozen PySCF NAC mapping from v0.23.2 is unchanged and is not exercised by this
fixed-geometry calculation.
