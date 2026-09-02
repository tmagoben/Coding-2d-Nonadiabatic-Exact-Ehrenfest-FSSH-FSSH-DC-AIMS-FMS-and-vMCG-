# v0.19 Detailed Derivations

## 1. Linear generalized-coordinate projection

Let flattened Cartesian coordinates be

$$
R
=
R_0+Jq.
$$

Then

$$
dR=J\,dq.
$$

For a scalar electronic energy,

$$
dE
=
(\nabla_R E)^T dR
=
(\nabla_R E)^T J\,dq.
$$

Therefore

$$
\boxed{
\nabla_q E
=
J^T\nabla_R E.
}
$$

The same covector transformation applies to derivative couplings:

$$
\boxed{
d_{IJ}^{(q)}
=
J^T d_{IJ}^{(R)}.
}
$$

For Cartesian mass matrix $M_R$,

$$
T
=
\frac12\dot R^T M_R\dot R
=
\frac12\dot q^T J^TM_RJ\dot q,
$$

so

$$
\boxed{
M_q=J^T M_R J.
}
$$

## 2. Cartesian lift used by the deterministic validation backend

The synthetic backend begins with known generalized covectors $g_q$ and must construct
a Cartesian covector $g_R$ satisfying

$$
J^Tg_R=g_q.
$$

For full-column-rank $J$, the minimum-subspace lift used by v0.19 is

$$
\boxed{
g_R
=
J(J^TJ)^{-1}g_q.
}
$$

Indeed,

$$
J^Tg_R
=
J^TJ(J^TJ)^{-1}g_q
=
g_q.
$$

The same construction is used for each NAC covector.

## 3. State-property transformation

Suppose tracked state $i$ corresponds to raw state $\pi(i)$ and

$$
|\phi_i'\rangle
=
p_i
|\phi_{\pi(i)}\rangle.
$$

Energies and diagonal gradients only reorder:

$$
E_i'=E_{\pi(i)},
$$

$$
\nabla E_i'
=
\nabla E_{\pi(i)}.
$$

Derivative couplings transform as

$$
\begin{aligned}
d_{ij}'
&=
\langle\phi_i'|\nabla\phi_j'\rangle
\\
&=
p_i^*p_j
\langle
\phi_{\pi(i)}
|
\nabla\phi_{\pi(j)}
\rangle
\end{aligned}
$$

for geometry-independent discrete phase corrections.

Thus

$$
\boxed{
d_{ij}'
=
p_i^*p_j
d_{\pi(i)\pi(j)}.
}
$$

## 4. Hungarian maximum-overlap assignment

Define

$$
W_{ij}=|O_{ij}|^2.
$$

The state-assignment problem is

$$
\boxed{
\pi_*
=
\arg\max_{\pi}
\sum_i W_{i,\pi(i)}.
}
$$

This is a maximum-weight perfect matching problem and is solved through a linear-sum
assignment algorithm.

The best assignment costs

$$
O(n_s^3).
$$

## 5. Exact second-best assignment without permutation enumeration

Let best assignment edges be

$$
B=
\{(i,\pi_*(i))\}.
$$

Any distinct perfect matching $A\neq B$ omits at least one edge of $B$.

For each best edge $e\in B$, solve the constrained problem

$$
S_e
=
\max_{A:e\notin A}S(A).
$$

Then every alternative matching belongs to at least one such constrained set.

Therefore

$$
\boxed{
S_{\rm second}
=
\max_{e\in B} S_e.
}
$$

This requires $n_s$ additional $O(n_s^3)$ assignments:

$$
\boxed{
O(n_s^4).
}
$$

The ambiguity margin remains

$$
\Delta S
=
S_{\rm best}-S_{\rm second}.
$$

## 6. Nearest-anchor tracking

Let the cache contain accepted generalized coordinates

$$
\mathcal Q
=
\{q_k\}.
$$

For a new point $q$,

$$
k_*
=
\arg\min_k
\|q-q_k\|_2.
$$

The raw new electronic snapshot is aligned against snapshot $k_*$.

This avoids a dependence on immediately preceding **call order**, but the current
implementation performs a linear nearest-cache search:

$$
\boxed{
O(N_{\rm cache}n_q)
}
$$

per cache miss.

Because an actual ab-initio electronic calculation is normally much more expensive
than this search for modest cache sizes, v0.19 keeps the transparent implementation.

A spatial index is a future large-cache optimization.

## 7. Root-label seed

At the first accepted geometry there is no previous electronic frame.

Therefore the first raw state order defines the tracked labels.

This is unavoidable without an externally specified electronic reference.

For reproducible branched dynamics, the intended workflow is

```text
evaluate initial physical geometry first
then allow center/centroid queries
```

The order-tolerance guarantee in v0.19 assumes that reference seed has been established.

## 8. Gauge graph overlap link

For two overlap-capable electronic snapshots,

$$
O_{uv}
=
\langle\Phi_u|\Phi_v\rangle.
$$

The discrete link is the nearest unitary polar factor

$$
O_{uv}
=
U\Sigma V^\dagger,
$$

$$
\boxed{
L_{uv}
=
UV^\dagger.
}
$$

Under local gauges

$$
\Phi_u\rightarrow\Phi_uG_u,
$$

$$
\Phi_v\rightarrow\Phi_vG_v,
$$

the link transforms covariantly:

$$
\boxed{
L_{uv}
\rightarrow
G_u^\dagger L_{uv}G_v.
}
$$

## 9. Pair-centroid reference

For pair $(i,j)$,

$$
q_c
=
\frac12(q_i+q_j).
$$

Both center electronic vectors are transported to $q_c$ before their scalar overlap
and Hamiltonian matrix element are evaluated.

This makes the pair approximation gauge covariant under consistent transformation of
the graph links and centroid operator.

It does not convert the discrete approximation into the complete continuous AIMS
kinetic-coupling integral.

## 10. Failure fallback

If a backend call fails at $q$ and nearest trusted cached point is $q_k$, fallback is
permitted only when

$$
\boxed{
\|q-q_k\|
\le
r_{\rm fallback}.
}
$$

The reused electronic quantities are marked stale in metadata.

The fallback is not inserted as a trusted gauge anchor.

Thus later tracking remains attached to actually evaluated electronic wavefunctions,
not synthetic fallback points.
