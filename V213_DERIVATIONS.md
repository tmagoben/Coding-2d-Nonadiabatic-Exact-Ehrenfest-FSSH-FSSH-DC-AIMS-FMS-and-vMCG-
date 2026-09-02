# v0.21.3 Detailed Derivations — SOC-Contract Freeze

v0.21.3 introduces no spin-orbit Hamiltonian. It freezes the mathematical meaning of
the objects a later SOC backend must supply and removes a gauge defect at exact
degeneracy.

## 1. Strict structural residuals

For matrices $A$ and $B$, define the explicit scaled Frobenius residual

$$
r(A,B)=
\frac{\lVert A-B\rVert_F}
{\max(\lVert A\rVert_F,\lVert B\rVert_F,1)}.
$$

The release checks

$$
r_H=r(H,H^\dagger),
\qquad
r_{K_a}=r(K_a,K_a^\dagger),
$$

$$
r_{D_a}=r(D_a,-D_a^\dagger),
\qquad
r_U=r(U^\dagger U,I).
$$

No implicit elementwise relative tolerance enters these identities. The default
operator tolerance is $10^{-12}$ and the default frame/isometry tolerance is
$10^{-10}$. The generalized mass matrix is also checked for explicit symmetry,
finiteness, and positive definiteness.

## 2. H, K, and D are distinct objects

For an electronic frame $\Phi(q)$, define

$$
H(q)=\Phi^\dagger\hat H_e(q)\Phi,
$$

$$
K_a(q)=\Phi^\dagger\left(\partial_a\hat H_e\right)\Phi,
$$

and

$$
D_a(q)=\Phi^\dagger\partial_a\Phi.
$$

The physical derivative operator $K_a$ is not generally the naive derivative of the
moving-frame matrix $H(q)$. Differentiating the matrix elements gives

$$
\boxed{
\partial_a H_{\rm matrix}=K_a+[H,D_a].
}
$$

Under $\Phi'=\Phi G(q)$ with $G\in U(s)$,

$$
H'=G^\dagger HG,
\qquad
K_a'=G^\dagger K_aG,
$$

$$
D_a'=G^\dagger D_aG+G^\dagger\partial_aG.
$$

Thus H and K transform covariantly, while D transforms as a connection. In a declared
fixed electronic frame, v0.21.3 requires $D_a=0$.

The future composition rule is already frozen:

$$
H=H_{\rm spin\text{-}free}+H_{\rm SOC},
$$

$$
\boxed{
K_a=K_{a,\rm spin\text{-}free}+K_{a,\rm SOC}.
}
$$

A backend cannot provide only an SOC energy matrix and silently omit its physical
derivative contribution. Nonzero SOC terms also require provenance with an explicit SOC
method. v0.21.3 exercises only the exactly zero-SOC branch.

## 3. Fixed electronic model space and provenance

The electronic model space declares, before propagation:

- a stable ordered state/component list;
- a fixed dimension $s$;
- fixed-frame or local-frame representation semantics;
- whether complete spin multiplets are required;
- hartree energy and bohr coordinate units;
- the rule that an electronic block is never pruned internally.

For complete multiplets, a root of multiplicity $m$ must contribute exactly $m$ unique
component labels. This prevents a later SOC calculation from silently changing the
meaning or dimension of a Gaussian's electronic block.

Canonical JSON containing the model, method versions, scalar-relativistic convention,
derivative method, SOC status, and numerical parameters is hashed with SHA-256:

$$
f=\operatorname{SHA256}(\operatorname{canonical\ provenance}).
$$

The fingerprint is metadata, a cache-key component, and a restart identity.

## 4. Degeneracy-safe density guidance

For a nonzero local coefficient vector $c_i$, define

$$
\rho_i=\frac{c_ic_i^\dagger}{c_i^\dagger c_i}.
$$

The generalized force policy is

$$
\boxed{
F_{i,a}=-\operatorname{Tr}(\rho_iK_a(q_i)).
}
$$

Under a local gauge,

$$
c_i'=G_i^\dagger c_i,
\qquad
\rho_i'=G_i^\dagger\rho_iG_i,
\qquad
K_a'=G_i^\dagger K_aG_i,
$$

so cyclicity of the trace gives $F_{i,a}'=F_{i,a}$.

When the local coefficient block falls below its amplitude threshold, choosing a lowest
eigenvector is ill-defined at degeneracy. Instead, let

$$
O_{01}=\Phi(q_0)^\dagger\Phi(q_1),
$$

and let $U_{01}$ be its nearest-unitary polar factor. Local coefficients transport as
$c_1\approx U_{01}^\dagger c_0$, so the guide density transports as

$$
\boxed{
\rho_1=U_{01}^\dagger\rho_0U_{01}.
}
$$

A populated block refreshes the guide density. A weak block retains the transported
density. A new zero block may inherit a parent's transported density. A genuinely
unseeded zero block receives zero force until a physical coefficient or explicit guide
density exists. No eigenvector is selected.

Predictor/corrector force calls are transactional: every trial starts from the last
accepted density checkpoint, rejected iterations are rolled back, and only the final
endpoint commits. Thus a near-threshold block cannot inherit an electronic density from
a numerically rejected trial.

## 5. Arbitrary-state fixed-frame projection

In a fixed global electronic frame,

$$
\Psi(q)=\sum_{i=1}^{N}\sum_{a=1}^{s}C_{ia}g_i(q)|a\rangle.
$$

The block metric is

$$
\boxed{
S=S^{\rm nuc}\otimes I_s,
}

where $S^{\rm nuc}_{ij}=\langle g_i|g_j\rangle$ is evaluated analytically for arbitrary
nuclear dimension and unequal positive-definite widths. Given a grid target
$\psi_a(q)$ and explicit quadrature weights $w(q)$,

$$
b_{ia}=\sum_q w(q)g_i(q)^*\psi_a(q),
$$

and the least-squares projection solves

$$
SC=b.
$$

No two-state or two-dimensional assumption remains. A local-frame calculation must
apply its explicit unitary frame transform; state labels are not inferred.

## 6. Provenance-safe complex cache

The fixed-frame cache key is a hash of

$$
(\text{format},\text{namespace},f,q),
$$

where $f$ is the complete provenance fingerprint. Each entry stores complex H, K, D,
state vectors, the real generalized mass matrix, and JSON metadata. Both the numerical
archive and metadata sidecar must exist; partial entries are rejected.

Cross-geometry overlaps in this cache must be the identity of the declared fixed frame.
Moving-frame many-electron wavefunction snapshots require a separate molecular restart
and overlap design and are refused rather than silently discarded.

## 7. Permanent spin-free branch

Setting

$$
H_{\rm SOC}=0,
\qquad
K_{a,\rm SOC}=0
$$

reproduces the spin-free H and K exactly. The Gaussian engine, block graph, propagation,
and convergence controls remain common to both future branches.
