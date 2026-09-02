# v0.24.2 PySCF connected-geometry SOC mathematics

This document fixes the exact mathematical path implemented in
`pyscf_differential_soc_v242.py`. Atomic units are used unless stated otherwise.

## 1. Common-orbital state average

At each geometry `q`, PySCF converges ROHF followed by an equal-weight
three-root SA-CASSCF(5e,4o) calculation. With root weights `w_I`,

$$
\sum_I w_I=1,
\qquad
D^{\mathrm{MO}}_{pq}(q)=\sum_I w_I D^{(I)}_{pq}(q).
$$

Doubly occupied inactive orbitals contribute two electrons and the active-space CI
roots contribute the remaining active density. The AO density is

$$
D^{\mathrm{AO}}(q)=C(q)D^{\mathrm{MO}}(q)C(q)^T,
$$

where v0.24.2 deliberately requires real common scalar orbitals and a real
spin-free density.

## 2. Direct-JK Breit-Pauli SOMF operator

For Cartesian SOC component `x`, `y`, or `z`, the one-electron AO tensor is the
nuclear-charge-weighted PySCF `int1e_prinvxp` integral. The two-electron mean field is

$$
G^{x}_{\mathrm{SOMF}}[D]
=J^x[D]-\frac32K^x_L[D]-\frac32K^x_R[D].
$$

PySCF's direct JK driver evaluates these three contractions from
`int2e_p1vxp1` shell blocks. The production path therefore never materializes the
rank-five array with shape `(3,n_AO,n_AO,n_AO,n_AO)`. The effective real AO tensor is

$$
h^{x}_{\mathrm{eff,AO}}
=\frac{1}{2c^2}\left(h^{x}_{1e}-G^{x}_{\mathrm{SOMF}}[D]\right),
$$

and the complex scalar-orbital tensor used by state interaction is

$$
h^{x}_{\mathrm{eff,MO}}
=-i\,C^T h^{x}_{\mathrm{eff,AO}}C.
$$

There is exactly one factor `1/(2c^2)`. At the central OH geometry, an explicit
rank-five contraction is evaluated only as a small-system oracle and agrees with the
direct-JK result to machine precision.

## 3. Spin-free roots and complete microstates

Each spin-free root `I` has energy `E_I`, spin `S_I`, and a reference CI vector at
`M_S=M_ref`. PySCF's `spin_square` verifies the declared doublet spin. Project-native
finite-sum Clebsch-Gordan algebra and Wigner-Eckart reconstruction form the direct SOC
matrix over the complete ordered states

$$
|I,S_I,M_S\rangle,
\qquad M_S=S_I,S_I-1,\ldots,-S_I.
$$

For the OH evidence this gives

$$
(D_1,+\tfrac12),(D_1,-\tfrac12),
(D_2,+\tfrac12),(D_2,-\tfrac12),
(D_3,+\tfrac12),(D_3,-\tfrac12).
$$

The state-interaction matrices are

$$
H(q)=H_{\mathrm{sf}}(q)+H_{\mathrm{SOC}}(q),
\qquad H_\alpha(q)=H_\alpha(q)^\dagger.
$$

## 4. Exact cross-geometry root overlaps

The existing restricted-CASSCF overlap engine embeds each active CI vector with its
doubly occupied core, evaluates determinant overlaps from the cross-geometry AO
overlap and orbital coefficients, and sums the determinant products to obtain

$$
O^{\mathrm{root}}_{IJ}(q_0,q)
=\langle\Psi_I(q_0)|\Psi_J(q)\rangle.
$$

This is a physical finite-manifold contraction, not generally a unitary matrix.
It is lifted conservatively into the complete microstate space:

$$
O_{I M,J M'}=
\delta_{S_I S_J}\,\delta_{M M'}\,O^{\mathrm{root}}_{IJ}.
$$

No overlap is invented between different spin or `M_S` sectors.

## 5. Degenerate-safe polar transport

For the singular-value decomposition

$$
O=U\Sigma V^\dagger,
$$

the right-to-left unitary transport is

$$
W=UV^\dagger.
$$

The raw contraction `O` remains evidence; only `W` transports coefficients and
operators. Singular-value retention, condition number, principal angle, contraction,
positive polar factor, and transport unitarity are independently audited. An endpoint
operator is represented in the center frame as

$$
\widetilde H_\alpha(q_\pm)
=W_\pm H_\alpha(q_\pm)W_\pm^\dagger.
$$

This full-subspace operation is invariant to arbitrary endpoint rotations inside the
near-degenerate first-two-root subspace. Independent root signs alone are not.

## 6. Component-resolved central differences

For OH bond displacement `h`,

$$
K_{\mathrm{sf}}^{(h)}
=\frac{\widetilde H_{\mathrm{sf}}(q_0+h)
-\widetilde H_{\mathrm{sf}}(q_0-h)}{2h},
$$

$$
K_{\mathrm{SOC}}^{(h)}
=\frac{\widetilde H_{\mathrm{SOC}}(q_0+h)
-\widetilde H_{\mathrm{SOC}}(q_0-h)}{2h},
$$

$$
K_{\mathrm{total}}^{(h)}
=K_{\mathrm{sf}}^{(h)}+K_{\mathrm{SOC}}^{(h)}.
$$

All four transported endpoint component matrices are serialized in each record, and
validation recomputes both component differences. Consequently a perturbation
`K_sf -> K_sf+Delta`, `K_SOC -> K_SOC-Delta` is rejected even though the stored total
is unchanged.

## 7. What the recorded connection means

The aligned contraction in the polar gauge is

$$
A_\pm=O_\pm W_\pm^\dagger.
$$

The code records

$$
D_{\mathrm{polar}}^{(h)}
=\operatorname{antiHerm}\left(\frac{A_+-A_-}{2h}\right),
\qquad
\operatorname{antiHerm}(X)=\frac{X-X^\dagger}{2}.
$$

The discarded Hermitian part is retained as a contraction/curvature diagnostic.
Because polar alignment chooses a local parallel-transport gauge, this construction
does not establish the continuous physical derivative connection of the SOC spinor
states. That capability claim remains false.

## 8. Convergence evidence

Centered differences have leading error `O(h^2)`, so on a truncation-dominated
plateau

$$
\frac{\|K^{(h/4)}-K^{(h/2)}\|_F}
{\|K^{(h/2)}-K^{(h)}\|_F}\longrightarrow\frac14.
$$

The frozen OH ladder is `h=(0.08,0.04,0.02)` bohr. Finer independent CASSCF solves
entered convergence noise and are not presented as superior evidence. The
Richardson estimate is

$$
K_R=\frac{4K^{(h/2)}-K^{(h)}}{3},
\qquad
\epsilon_R\approx\frac{\|K^{(h/2)}-K^{(h)}\|_F}{3}.
$$

This is a one-coordinate implementation/convergence preview, not a full molecular
derivative or accuracy benchmark.
