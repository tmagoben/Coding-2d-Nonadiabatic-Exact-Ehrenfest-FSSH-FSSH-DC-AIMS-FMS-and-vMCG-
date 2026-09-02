# v0.8 Detailed Derivations

## A. Moving-basis norm conservation

From

$$
iS\dot C=(H-iT)C,
$$

obtain

$$
S\dot C=-iHC-TC.
$$

Hermitian conjugation gives

$$
\dot C^\dagger S=iC^\dagger H-C^\dagger T^\dagger.
$$

Then

$$
\frac{d}{dt}(C^\dagger SC)
=C^\dagger[-T^\dagger+\dot S-T]C.
$$

Thus

$$
\boxed{
\dot S=T+T^\dagger
}
$$

is sufficient for exact norm conservation.

---

## B. Discrete metric correction

Let

$$
D_S=\frac{S_{n+1}-S_n}{\Delta t}.
$$

Suppose a physically motivated seed $T_0$ is known from the nuclear Gaussian motion.

Define

$$
\Delta T
=\frac12(D_S-T_0-T_0^\dagger).
$$

Since $D_S$ is Hermitian,

$$
\Delta T^\dagger=\Delta T.
$$

Set

$$
T=T_0+\Delta T.
$$

Then

$$
T+T^\dagger
=T_0+T_0^\dagger+2\Delta T
=D_S.
$$

Therefore

$$
\boxed{
\frac{S_{n+1}-S_n}{\Delta t}=T+T^\dagger.
}
$$

Because $\Delta T$ is Hermitian, the anti-Hermitian part of $T_0$ is unchanged.

---

## C. Temporal overlap and derivative coupling

Expand

$$
|\phi_j(t+\Delta t)\rangle
=|\phi_j(t)\rangle+\Delta t|\dot\phi_j(t)\rangle+\mathcal O(\Delta t^2).
$$

Then

$$
O_{ij}
=\langle\phi_i(t)|\phi_j(t+\Delta t)\rangle
=\delta_{ij}+\Delta t\langle\phi_i|\dot\phi_j\rangle+\mathcal O(\Delta t^2).
$$

Along $R(t)$,

$$
|\dot\phi_j\rangle
=\dot R\cdot\nabla_R|\phi_j\rangle,
$$

so

$$
\boxed{
O=I+\Delta t(\dot R\cdot d)+\mathcal O(\Delta t^2).
}
$$

For an orthonormal basis $d$ is anti-Hermitian, so the first-order overlap change is itself an infinitesimal unitary rotation.

---

## D. Local-diabatic overlap step

Start from

$$
c_{n+1}
=e^{-iH_{n+1}\Delta t/2}L^\dagger e^{-iH_n\Delta t/2}c_n.
$$

Expand each factor:

$$
e^{-iH_n\Delta t/2}=I-\frac{i\Delta t}{2}H_n+\mathcal O(\Delta t^2),
$$

$$
L^\dagger=I-\Delta t(\dot R\cdot d)+\mathcal O(\Delta t^2),
$$

$$
e^{-iH_{n+1}\Delta t/2}=I-\frac{i\Delta t}{2}H_{n+1}+\mathcal O(\Delta t^2).
$$

For a smooth Hamiltonian,

$$
H_{n+1}=H_n+\mathcal O(\Delta t).
$$

Multiplying to first order,

$$
c_{n+1}
=\left[I-iH_n\Delta t-(\dot R\cdot d)\Delta t\right]c_n+\mathcal O(\Delta t^2).
$$

Thus

$$
\boxed{
\dot c=-iHc-(\dot R\cdot d)c.
}
$$

---

## E. Gauge covariance of the polar link

If

$$
O'=G_n^\dagger O G_{n+1},
$$

with unitary endpoint gauges, and

$$
O=U\Sigma V^\dagger,
$$

then one SVD of $O'$ is

$$
O'=(G_n^\dagger U)\Sigma(G_{n+1}^\dagger V)^\dagger.
$$

Its polar factor is

$$
L'
=(G_n^\dagger U)(G_{n+1}^\dagger V)^\dagger
=G_n^\dagger UV^\dagger G_{n+1}.
$$

Therefore

$$
\boxed{
L'=G_n^\dagger L G_{n+1}.
}
$$

---

## F. Generalized-mass spawn momentum

Let

$$
p_b=p_a+\lambda n,
$$

with inverse mass matrix $B=M^{-1}$.

Energy conservation is

$$
\frac12(p_a+\lambda n)^TB(p_a+\lambda n)+E_b
=\frac12p_a^TBp_a+E_a.
$$

Cancel the common term:

$$
(n^TBn)\lambda^2+2(p_a^TBn)\lambda+2(E_b-E_a)=0.
$$

Hence

$$
\boxed{
\lambda
=\frac{-p_a^TBn\pm\sqrt{(p_a^TBn)^2-2(n^TBn)(E_b-E_a)}}{n^TBn}.
}
$$

The implementation chooses the real root with the smaller magnitude.

---

## G. Basis insertion continuity

Before spawning,

$$
\Psi^- = \sum_{i=1}^{N} C_iG_i.
$$

After inserting $G_{N+1}$ with coefficient zero,

$$
\Psi^+=\sum_{i=1}^{N}C_iG_i+0\,G_{N+1}.
$$

Therefore

$$
\boxed{
\Psi^+=\Psi^-.
}
$$

The physical state changes only when later propagation makes $C_{N+1}\neq0$.

---

## H. Graph cycle rank

For a connected finite graph, choose a spanning tree. The tree contains

$$
N_V-1
$$

edges.

Every additional edge creates one independent fundamental cycle. Therefore

$$
\boxed{
N_{\mathrm{cycle}}=N_E-(N_V-1)=N_E-N_V+1.
}
$$

This number increases naturally as TBF-pair centroid paths are added through time.
