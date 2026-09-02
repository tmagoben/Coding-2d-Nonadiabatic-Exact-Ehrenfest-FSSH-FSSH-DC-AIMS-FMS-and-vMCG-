# v0.23.2 PySCF NAC convention erratum

## Internal convention

For electronic states indexed by i and j, the dynamics convention is

$$
d_{ij}(R)=\langle\Phi_i(R)|\nabla_R\Phi_j(R)\rangle,
\qquad d_{ji}=-d_{ij}^{*}.
$$

The physical off-diagonal Hamiltonian derivative satisfies

$$
(K_a)_{ij}=(E_j-E_i)(d_a)_{ij}.
$$

## Correction

PySCF documentation describes `state=(ket,bra)` as returning
`<bra|d ket/dR>`. A literal reading led v0.23.1 and earlier adapters to request
`state=(j,i)` for internal `d[i,j]`. Direct PySCF 2.13.1 execution, tested against
phase-aligned many-electron overlap central differences, shows that the
production mapping required by this framework is instead:

```text
internal d[i,j]  <-  PySCF state=(i,j), mult_ediff=False, use_etfs=False
```

Requesting `(j,i)` gives the antisymmetric counterpart and therefore the
opposite sign for a real pair. All production PySCF paths now use the centralized
`pyscf_state_tuple_for_internal_dij_v232(i, j)` helper.

Root phases may change with the eigensolver or BLAS runtime. The certification
therefore aligns displaced-state phases before differentiating overlaps and does
not hard-code an absolute, unaligned root sign.

## ETF distinction

`use_etfs=False` matches the derivative of the full many-electron overlap.
`use_etfs=True` removes the translational component and is scientifically useful,
but it is not interchangeable with the full-overlap derivative used by the
framework contract.

This document supersedes the adapter mapping stated in `V03_PYSCF_NOTES.md`,
`V05_PYSCF_BACKEND.md`, and the historical v0.23.1 adapter notes.
