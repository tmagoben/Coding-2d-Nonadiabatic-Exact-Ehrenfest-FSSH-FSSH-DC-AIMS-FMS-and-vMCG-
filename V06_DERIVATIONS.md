# v0.6 Detailed Derivations

---

## A. Nonorthogonal Slater determinant overlap

Let

$$
|\mathcal D_A\rangle
=
a_1^\dagger\cdots a_N^\dagger|0\rangle
$$

and

$$
|\mathcal D_B\rangle
=
b_1^\dagger\cdots b_N^\dagger|0\rangle.
$$

Define the one-particle overlap matrix

$$
S_{ij}
=
\langle a_i|b_j\rangle.
$$

Using fermionic antisymmetry, the $N$-electron overlap is

$$
\boxed{
\langle\mathcal D_A|\mathcal D_B\rangle
=
\det S.
}
$$

For separated alpha and beta determinants,

$$
|\mathcal D\rangle
=
|\mathcal D_\alpha\rangle
|\mathcal D_\beta\rangle,
$$

the spin sectors factor:

$$
\boxed{
\langle\mathcal D_A|\mathcal D_B\rangle
=
\det S_\alpha
\det S_\beta.
}
$$

---

## B. Multiconfigurational overlap

Let

$$
|\Psi_A\rangle
=
\sum_I C_I^A|D_I^A\rangle,
$$

$$
|\Psi_B\rangle
=
\sum_J C_J^B|D_J^B\rangle.
$$

Then

$$
\langle\Psi_A|\Psi_B\rangle
=
\sum_{IJ}
(C_I^A)^*
C_J^B
\langle D_I^A|D_J^B\rangle.
$$

Therefore

$$
\boxed{
\langle\Psi_A|\Psi_B\rangle
=
\sum_{IJ}
(C_I^A)^*
C_J^B
\det S_{IJ}.
}
$$

For spin-separated CI coefficients $C_{KL}$,

$$
\boxed{
\langle\Psi_A|\Psi_B\rangle
=
\sum_{KL,MN}
(C_{KL}^A)^*
C_{MN}^B
\det S_{\alpha;KM}
\det S_{\beta;LN}.
}
$$

This is what the PySCF FCI overlap helper evaluates after the CASSCF CI vector is
embedded into the core+active determinant space.

---

## C. Core+active determinant embedding

Suppose the active alpha determinant is encoded by bit string

$$
b_{\mathrm{act}}.
$$

Let the $n_c$ core orbitals be the lowest orbital indices.

Their occupied bit mask is

$$
\boxed{
b_{\mathrm{core}}=2^{n_c}-1.
}
$$

Shift the active string by $n_c$ positions:

$$
b_{\mathrm{act}}^{\mathrm{shifted}}
=
b_{\mathrm{act}}\ll n_c.
$$

Then the full core+active determinant string is

$$
\boxed{
b_{\mathrm{full}}
=
b_{\mathrm{core}}
\lor
(b_{\mathrm{act}}\ll n_c).
}
$$

The same construction is performed independently for alpha and beta strings.

The CI amplitude is copied into the corresponding determinant address of the
core+active FCI vector.

All other full-space determinants retain zero coefficient.

---

## D. Example showing why core-active cross blocks matter

Take one core spatial orbital and one active spatial orbital.

Assume the current correlated orbitals are a rotation of the previous ones:

$$
S_{\mathrm{MO}}
=
\begin{pmatrix}
\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta
\end{pmatrix}.
$$

Let the active electron count be one alpha electron and zero beta electrons.

The alpha determinant occupies both core and active orbitals, so

$$
S_\alpha
=
S_{\mathrm{MO}}
$$

and

$$
\det S_\alpha=1.
$$

The beta determinant occupies only the core orbital, so

$$
\det S_\beta=\cos\theta.
$$

Therefore the exact many-electron overlap is

$$
\boxed{
\langle\Psi_A|\Psi_B\rangle
=
\cos\theta.
}
$$

A naive separated core/active product would generally give a different result because
it discards the off-diagonal core-active overlap.

This exact case appears in the v0.6 test suite.

---

## E. Maximum-overlap assignment

For permutation $\pi$, define

$$
\mathcal S(\pi)
=
\sum_i|O_{i,\pi(i)}|^2.
$$

The tracker chooses

$$
\boxed{
\pi^*
=
\arg\max_\pi\mathcal S(\pi).
}
$$

This is a linear assignment problem.

Because the intended state manifolds are small, v0.6 enumerates all permutations.
For large $N$, one would replace this with a Hungarian/linear-sum assignment solver.

---

## F. Gauge correction

Let

$$
z_i=O_{i,\pi(i)}.
$$

Choose current ket phase $p_i$.

The transformed overlap is

$$
\langle\Psi_i^{\mathrm{prev}}|p_i\Psi_{\pi(i)}^{\mathrm{curr}}\rangle
=
p_i z_i.
$$

To make this positive real,

$$
p_i z_i=|z_i|,
$$

so

$$
\boxed{
p_i=\frac{z_i^*}{|z_i|}.
}
$$

For real states,

$$
\boxed{
p_i=\operatorname{sign}(z_i).
}
$$

---

## G. NAC gauge transformation

Start from

$$
d_{ij}
=
\langle\Phi_i|\nabla\Phi_j\rangle.
$$

Transform

$$
|\Phi_i'\rangle=p_i|\Phi_i\rangle
$$

with geometry-independent phase at that discrete tracking step.

Then

$$
\langle\Phi_i'|
=
p_i^*\langle\Phi_i|.
$$

Therefore

$$
d'_{ij}
=
p_i^*p_j
\langle\Phi_i|\nabla\Phi_j\rangle.
$$

Hence

$$
\boxed{
d'_{ij}=p_i^*p_jd_{ij}.
}
$$

For $p_i=\pm1$,

$$
\boxed{
d'_{ij}=p_ip_jd_{ij}.
}
$$

---

## H. Procrustes alignment

Let

$$
O=\Phi_A^\dagger\Phi_B.
$$

We seek unitary $Q$ minimizing

$$
\|\Phi_BQ-\Phi_A\|_F^2.
$$

Expanding and using orthonormality gives the equivalent maximization

$$
\max_Q
\operatorname{Re}
\operatorname{Tr}(OQ).
$$

Take

$$
O=U\Sigma V^\dagger.
$$

The optimum is

$$
\boxed{
Q=VU^\dagger.
}
$$

Then

$$
OQ
=
U\Sigma U^\dagger,
$$

which is Hermitian positive semidefinite.

---

## I. Directional NAC from overlap

Expand

$$
|\Phi_j(s+\Delta s)\rangle
=
|\Phi_j(s)\rangle
+
\Delta s|\partial_s\Phi_j\rangle
+
\frac12\Delta s^2|\partial_s^2\Phi_j\rangle+\cdots.
$$

Project with $\langle\Phi_i(s)|$:

$$
O_{ij}
=
\delta_{ij}
+
\Delta s d_{ij}
+
\mathcal O(\Delta s^2).
$$

Take the Hermitian conjugate:

$$
O^\dagger
=
I
+
\Delta s d^\dagger
+
\mathcal O(\Delta s^2).
$$

Because

$$
d^\dagger=-d,
$$

$$
O-O^\dagger
=
2\Delta s\,d
+
\mathcal O(\Delta s^2).
$$

Therefore

$$
\boxed{
d
=
\frac{O-O^\dagger}{2\Delta s}
+
\mathcal O(\Delta s).
}
$$

---

## J. Principal angles

For subspace overlap matrix

$$
O_{\mathrm{sub}}
=
\Phi_A^\dagger\Phi_B,
$$

let its singular values be

$$
\sigma_k.
$$

The principal angles are

$$
\boxed{
\vartheta_k
=
\arccos\sigma_k.
}
$$

If all $\sigma_k\approx1$, the two subspaces are nearly the same even if their
individual basis vectors have rotated strongly.

This is the correct diagnostic near an exact or near degeneracy.

---

## K. Topological obstruction around a CI

Suppose a real electronic eigenvector is parallel transported continuously around a
closed path enclosing one conical intersection.

Locally one can always choose adjacent signs so

$$
\langle\Phi(s)|\Phi(s+\Delta s)\rangle>0.
$$

But after a complete loop,

$$
\boxed{
|\Phi_{\mathrm{final}}\rangle
=
-|\Phi_{\mathrm{initial}}\rangle.
}
$$

Thus there is no globally single-valued real phase convention over a domain enclosing
the degeneracy.

A state tracker should preserve this holonomy rather than force the final sign back to
$+1$ and erase the geometric phase.
