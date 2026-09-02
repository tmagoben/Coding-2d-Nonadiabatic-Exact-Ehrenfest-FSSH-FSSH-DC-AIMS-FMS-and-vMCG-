# v0.21.4 Detailed Derivations

v0.21.4 certifies that an electronic provider describes one differentiable model and
that a self-consistent trajectory can resume without losing hidden algorithmic state.
It introduces no physical spin-orbit term.

## 1. Why pointwise validation is insufficient

At one geometry, the conditions

$$
H=H^\dagger,\qquad K_a=K_a^\dagger,\qquad D_a=-D_a^\dagger
$$

do not prove that K and D are derivatives of the supplied electronic model. A provider
can pass all three identities while returning an arbitrary Hermitian K or an arbitrary
anti-Hermitian D. Certification must therefore compare neighboring geometries.

Let

$$
O_{0\pm}=\langle\Phi(q)|\Phi(q\pm h_a e_a)\rangle,
$$

and let $U_{0\pm}$ be the nearest-unitary polar factors of these overlaps. Neighboring
Hamiltonians transported to the center frame are

$$
H_{\pm\to0}=U_{0\pm}H(q\pm h_a e_a)U_{0\pm}^\dagger.
$$

The centered physical-operator derivative estimate is

$$
\boxed{
K_a^{\mathrm{FD}}(q)=\frac{H_{+\to0}-H_{-\to0}}{2h_a}.
}
$$

Transport removes the frame derivative before differencing. Consequently this object
is compared with the physical derivative operator K, rather than with the derivative of
the raw moving-frame matrix H.

The derivative connection is independently estimated from the raw overlaps:

$$
\boxed{
D_a^{\mathrm{FD}}(q)=\frac{O_{0+}-O_{0-}}{2h_a}.
}
$$

For either comparison the scaled Frobenius error is

$$
\epsilon(A,B)=
\frac{\lVert A-B\rVert_F}
{\max(\lVert A\rVert_F,\lVert B\rVert_F,1)}.
$$

The audit also checks overlap isometry, pointwise H/K/D/mass structure, exact displaced
geometry identity, finite data, and one unchanged provenance fingerprint.

## 2. Truncation and step selection

Centered differences have $O(h_a^2)$ truncation error, while floating-point
cancellation grows as $h_a$ becomes too small. The default $h_a=10^{-4}$ is a
deterministic analytic-fixture value, not a universal molecular-backend choice. A real
backend must document coordinate units, numerical noise, and a step-refinement plateau.
Per-coordinate positive steps are therefore part of the audit settings.

## 3. Exact zero-SOC rehearsal

The future operator decomposition is

$$
H=H_{\mathrm{sf}}+H_{\mathrm{SOC}},\qquad
K_a=K_{a,{\mathrm{sf}}}+K_{a,{\mathrm{SOC}}}.
$$

v0.21.4 sets

$$
H_{\mathrm{SOC}}=0,\qquad K_{a,{\mathrm{SOC}}}=0
$$

as explicit complex arrays and routes them through the same composition function a
later optional physical backend will use. Exact equality is required for H, K, D, the
mass matrix, and cross-geometry overlaps. Dynamics is then compared in the Gaussian
metric after global-phase alignment. This rehearses integration plumbing only; it is
not a physical SOC calculation.

## 4. Complete restart state

For N Gaussian trajectory basis functions, electronic dimension s, and nuclear
dimension d, a v0.21.4 checkpoint stores

$$
\{u_i,q_i,p_i,A_i\}_{i=1}^{N},\qquad C\in\mathbb{C}^{Ns},
$$

together with accepted density-guide state

$$
\{m_i,\rho_i\}_{i=1}^{N},
$$

six guidance counters, and the canonical active sparse-edge set

$$
\mathcal E_{\mathrm{UID}}=\{(\min(u_i,u_j),\max(u_i,u_j))\}.
$$

Stable UIDs are essential: array indices can change after insertion or pruning, whereas
sparse hysteresis belongs to physical Gaussian identities.

The manifest stores format, step, time, dt, s, d, provider fingerprint, and settings
fingerprint. Time must satisfy

$$
t=n\,\Delta t
$$

within an explicit scaled tolerance.

## 5. Canonical identity and integrity

Propagation settings are serialized canonically together with the release semantics and
hashed:

$$
f_{\mathrm{settings}}=\operatorname{SHA256}(\operatorname{canonical\ settings}).
$$

The checkpoint integrity digest consumes the canonical manifest, then for each array in
a fixed order consumes its name, exact dtype, shape, and contiguous bytes:

$$
f_{\mathrm{checkpoint}}=\operatorname{SHA256}
\left(M\,\Vert\,\big\Vert_k(n_k,\tau_k,s_k,b_k)\right).
$$

This detects both metadata and numerical changes. It is an integrity mechanism, not an
authentication or secrecy mechanism.

## 6. Deterministic continuation

At resume, the checkpoint is validated against the declared provider and settings.
Density guides and sparse UID edges are restored before the first continued step.
Provider snapshots are rebuilt at the stored geometries, avoiding unsafe serialization
of backend-specific wavefunction objects.

If the checkpoint ended at step $n_0$, a local adaptation callback at segment step k is
called with global step

$$
n=n_0+k.
$$

Stored records, corrector histories, and adaptation events receive the same global
offset. The release compares uninterrupted and segmented final positions, momenta, and
metric-aware phase-aligned coefficients.

