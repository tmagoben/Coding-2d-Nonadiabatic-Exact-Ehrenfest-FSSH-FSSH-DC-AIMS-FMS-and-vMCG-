# v0.8 Theory: Time-Dependent Graph Gaussian Dynamics

Version 0.8 turns the static gauge graph of v0.7 into a **time-dependent graph of moving Gaussian trajectory basis functions (TBFs)**.

The central problem is now no longer only

$$
\text{How do we compare electronic states at different geometries?}
$$

but

$$
\boxed{
\text{How do we propagate a nonorthogonal Gaussian wavefunction while its electronic gauge graph itself grows in time?}
}
$$

The release therefore combines:

1. temporal electronic-overlap links;
2. explicit-NAC and overlap/local-diabatic electronic propagation;
3. moving nonorthogonal Gaussian matrices $S(t)$, $H(t)$, and $T(t)$;
4. dynamic zero-amplitude spawning;
5. incremental TBF-center and pair-centroid graph growth;
6. a PySCF many-electron-snapshot graph interface.

Atomic units are used.

---

## 1. Time-dependent Gaussian expansion

Write

$$
\boxed{
|\Psi(t)\rangle
=
\sum_{k=1}^{N_G(t)} C_k(t)|G_k(t)\rangle.
}
$$

The number of TBFs may change after a spawning event.

Each basis function contains a nuclear Gaussian and a local electronic state:

$$
|G_k(t)\rangle
=
|g_k(\mathbf q,t)\rangle
|e_k(\mathbf R_k,t)\rangle.
$$

Because both factors move,

$$
\frac{d}{dt}|G_k\rangle
=
|\dot g_k\rangle|e_k\rangle
+
|g_k\rangle|\dot e_k\rangle.
$$

The second term is precisely where electronic gauge continuity enters the moving basis.

---

## 2. Projected TDSE for a moving basis

Define

$$
S_{ij}=\langle G_i|G_j\rangle,
$$

$$
H_{ij}=\langle G_i|\hat H|G_j\rangle,
$$

and

$$
\boxed{
T_{ij}=\langle G_i|\dot G_j\rangle.
}
$$

Substitute the expansion into

$$
i|\dot\Psi\rangle=\hat H|\Psi\rangle.
$$

After projection with $\langle G_i|$,

$$
\boxed{
iS\dot C=(H-iT)C.
}
$$

This is the fundamental coefficient equation propagated in v0.8.

---

## 3. Metric compatibility

Differentiate the overlap matrix:

$$
\dot S_{ij}
=
\langle\dot G_i|G_j\rangle
+
\langle G_i|\dot G_j\rangle.
$$

Since

$$
\langle\dot G_i|G_j\rangle=T_{ji}^*,
$$

we obtain

$$
\boxed{
\dot S=T+T^\dagger.
}
$$

This identity is essential. It is not merely a numerical convenience.

---

## 4. Norm conservation in a moving nonorthogonal basis

The physical norm is

$$
\boxed{
N=C^\dagger S C.
}
$$

Differentiate:

$$
\dot N
=
\dot C^\dagger SC
+C^\dagger\dot S C
+C^\dagger S\dot C.
$$

From

$$
iS\dot C=(H-iT)C,
$$

we have

$$
S\dot C=-iHC-TC.
$$

The Hermitian conjugate is

$$
\dot C^\dagger S
=iC^\dagger H-C^\dagger T^\dagger.
$$

Therefore

$$
\dot N
=
C^\dagger
[-T^\dagger+\dot S-T]
C.
$$

Using

$$
\dot S=T+T^\dagger,
$$

we obtain

$$
\boxed{
\dot N=0.
}
$$

Thus a correct moving-basis propagation must treat $S$ and $T$ consistently.

---

# 5. Electronic propagation in an adiabatic representation

Along a nuclear path $\mathbf R(t)$,

$$
|\psi_e\rangle
=
\sum_a c_a(t)|\phi_a(\mathbf R(t))\rangle.
$$

The adiabatic coefficient equation is

$$
\boxed{
\dot c_a
=
-iE_a c_a
-
\sum_b
\dot{\mathbf R}\cdot\mathbf d_{ab}\,c_b.
}
$$

Define

$$
D_{ab}
=
\dot{\mathbf R}\cdot\mathbf d_{ab}.
$$

For a real adiabatic basis,

$$
D^T=-D.
$$

Hence

$$
\boxed{
H_{\rm eff}=E-iD
}
$$

is Hermitian.

v0.8 propagates the explicit-NAC reference through

$$
\boxed{
c(t+\Delta t)
=
\exp[-iH_{\rm eff}(t+\Delta t/2)\Delta t]c(t).
}
$$

---

# 6. Temporal electronic overlap

Now consider two successive electronic frames,

$$
\Phi_n
=
\{\phi_a(\mathbf R_n)\},
$$

$$
\Phi_{n+1}
=
\{\phi_a(\mathbf R_{n+1})\}.
$$

Their overlap matrix is

$$
\boxed{
O_{n,n+1}
=
\Phi_n^\dagger\Phi_{n+1}.
}
$$

The raw overlap need not be exactly unitary for a finite selected electronic subspace. v0.8 uses its unitary polar factor,

$$
O=U\Sigma V^\dagger,
$$

$$
\boxed{
L=UV^\dagger.
}
$$

$L$ is the discrete electronic connection carried by the temporal graph edge.

---

# 7. Coefficient transport between electronic frames

If

$$
|\psi\rangle
=\Phi_n c_n
=\Phi_{n+1}c_{n+1},
$$

and

$$
L\approx\Phi_n^\dagger\Phi_{n+1},
$$

then

$$
\boxed{
c_{n+1}\approx L^\dagger c_n.
}
$$

This transformation is not a stochastic hop. It is simply the representation change between neighboring electronic frames.

---

# 8. Overlap/local-diabatic Strang propagator

v0.8 combines basis transport and electronic phase evolution symmetrically:

$$
\boxed{
c_{n+1}
=
e^{-iH_{n+1}\Delta t/2}
L_{n,n+1}^\dagger
e^{-iH_n\Delta t/2}
c_n.
}
$$

This is implemented in `temporal_electronic.py`.

It has three useful properties:

1. each factor is unitary;
2. the whole step is explicitly gauge covariant;
3. it converges to the explicit derivative-coupling equation for small steps.

---

# 9. Small-step relation to the NAC equation

For a small displacement,

$$
\mathbf R_{n+1}
=
\mathbf R_n
+
\dot{\mathbf R}\Delta t,
$$

expand

$$
|\phi_b(\mathbf R_{n+1})\rangle
=
|\phi_b(\mathbf R_n)\rangle
+
\Delta t\,
\dot{\mathbf R}\cdot\nabla|\phi_b\rangle
+
\mathcal O(\Delta t^2).
$$

Therefore

$$
O_{ab}
=
\delta_{ab}
+
\Delta t\,
\dot{\mathbf R}\cdot\mathbf d_{ab}
+
\mathcal O(\Delta t^2).
$$

Hence

$$
L^\dagger
=
I
-
\Delta t\,
\dot{\mathbf R}\cdot d
+
\mathcal O(\Delta t^2).
$$

Together with

$$
e^{-iH\Delta t}
=
I-iH\Delta t+\mathcal O(\Delta t^2),
$$

the overlap step yields

$$
c_{n+1}
=
\left[
I
-iE\Delta t
-
(\dot{\mathbf R}\cdot d)\Delta t
\right]c_n
+
\mathcal O(\Delta t^2),
$$

which is the adiabatic electronic equation.

---

# 10. Local gauge covariance of the overlap step

Let the two endpoint frames be independently transformed:

$$
\Phi_n' = \Phi_n G_n,
$$

$$
\Phi_{n+1}' = \Phi_{n+1}G_{n+1}.
$$

Then

$$
H_n'=G_n^\dagger H_nG_n,
$$

$$
H_{n+1}'=G_{n+1}^\dagger H_{n+1}G_{n+1},
$$

and

$$
O' = G_n^\dagger O G_{n+1}.
$$

The polar link transforms covariantly:

$$
\boxed{
L'=G_n^\dagger L G_{n+1}.
}
$$

Since

$$
c_n'=G_n^\dagger c_n,
$$

the propagated state satisfies

$$
\boxed{
c_{n+1}'=G_{n+1}^\dagger c_{n+1}.
}
$$

The v0.8 test suite checks this numerically using independent random $U(2)$ gauges at the two time slices.

---

# 11. Time-dependent electronic graph

In v0.7, graph nodes were static geometries.

In v0.8 each active TBF acquires a new electronic-frame node after every nuclear step:

```text
TBF 0:  q0(t0) -- q0(t1) -- q0(t2) -- ...
TBF 1:             q1(t1) -- q1(t2) -- ...
```

The temporal edges store electronic overlap links.

For every interacting TBF pair, a centroid node is also added:

```text
q_i(t) ---- centroid_ij(t) ---- q_j(t)
```

so the Hamiltonian matrix element uses one common local electronic reference frame.

The graph therefore grows in both **time** and **branching structure**.

---

# 12. Graph cycle count

For a connected graph with

$$
N_V
$$

vertices and

$$
N_E
$$

edges, the cycle rank is

$$
\boxed{
\beta_1=N_E-N_V+1.
}
$$

Every new center-centroid-center path can create additional loops.

Those loops are not redundant bookkeeping: their Wilson products diagnose whether the discrete electronic connection is globally flat or carries nontrivial holonomy.

---

# 13. Moving Gaussian nuclear contribution to $T$

For a fixed-width Gaussian,

$$
g_j
=
N\exp
\left[
-\frac12(q-q_j)^TA(q-q_j)
+ip_j^T(q-q_j)
\right],
$$

we already derived

$$
\dot g_j
=
\left[
(A(q-q_j)-ip_j)^T\dot q_j
+i(q-q_j)^T\dot p_j
\right]g_j.
$$

Hence the exact nuclear moving-basis matrix element is

$$
\boxed{
T_{ij}^{\rm nuc}
=\langle g_i|\dot g_j\rangle.
}
$$

In a common graph-transported electronic frame, v0.8 multiplies this by the transported electronic overlap.

This supplies the **seed connection**.

---

# 14. Discrete metric-compatible electronic correction

Finite graph updates determine

$$
S_n
$$

and

$$
S_{n+1}.
$$

Approximate

$$
\dot S
\approx
\frac{S_{n+1}-S_n}{\Delta t}.
$$

Let

$$
T^{(0)}
$$

be the exact nuclear seed.

Its Hermitian part may not contain the full electronic-frame motion. v0.8 therefore defines

$$
\boxed{
T
=
T^{(0)}
+
\frac12
\left[
\dot S
-T^{(0)}
-T^{(0)\dagger}
\right].
}
$$

Then

$$
T+T^\dagger
=
\dot S
$$

exactly at the discrete level.

The correction is Hermitian, so the anti-Hermitian part of the physically derived nuclear seed is preserved.

This is the **minimal metric-compatible electronic connection approximation** used in the v0.8 prototype.

---

# 15. Midpoint coefficient propagation

The old and new matrices are

$$
S_n,H_n,
$$

and

$$
S_{n+1},H_{n+1}.
$$

Use

$$
S_{1/2}=\frac12(S_n+S_{n+1}),
$$

$$
H_{1/2}=\frac12(H_n+H_{n+1}).
$$

The local coefficient generator is

$$
\boxed{
A_C
=
S_{1/2}^{-1}
(-iH_{1/2}-T_{1/2}).
}
$$

v0.8 propagates

$$
\dot C=A_CC
$$

with RK4 over one small interval.

No explicit $S^{-1}$ matrix is formed in the coefficient equation itself; the generator is obtained by a linear solve.

---

# 16. Dynamic spawning

For a parent TBF on state $a$, the coupling indicator remains

$$
\boxed{
\eta_{ab}
=
|\dot q^Td_{ab}|.
}
$$

If

$$
\eta_{ab}>\eta_{\rm spawn},
$$

a child candidate on state $b$ is constructed.

The momentum update along the NAC direction is the generalized-mass energy-conserving solution already derived in v0.5.

---

# 17. Spawn insertion occurs between time intervals

A spawn is inserted only after the current old-basis propagation interval is finished.

Before insertion,

$$
\Psi
=
\sum_{i=1}^NC_iG_i.
$$

After adding the child,

$$
\boxed{
C'=(C_1,\ldots,C_N,0)^T.
}
$$

Therefore

$$
\boxed{
\Psi_{\rm after}=\Psi_{\rm before}.
}
$$

The new amplitude becomes nonzero only during later coupled propagation.

---

# 18. Child and parent electronic frame at birth

At the instant of spawning, parent and child occupy the same nuclear position.

Therefore they can share the same **electronic graph node** while carrying different electronic coefficient vectors:

$$
e_a=(0,\ldots,1_a,\ldots)^T,
$$

$$
e_b=(0,\ldots,1_b,\ldots)^T.
$$

On the next nuclear step their phase-space centers diverge, after which each receives a separate time-labelled graph node.

This avoids inventing two distinct electronic calculations at an identical geometry.

---

# 19. Explicit NAC versus overlap propagation benchmark

The v0.8 CI benchmark propagates the same electronic state along the same path using:

### Method A

$$
\boxed{
H_{\rm eff}=E-i\dot R\cdot d.
}
$$

### Method B

$$
\boxed{
e^{-iH_{n+1}\Delta t/2}
L_{n,n+1}^\dagger
e^{-iH_n\Delta t/2}.
}
$$

At 400 steps for the supplied path, the phase-insensitive fidelity is approximately

$$
\boxed{
F\approx0.9999999999996.
}
$$

Thus the two implementations agree to essentially machine-visible precision on the refined analytic benchmark.

This is one of the most important cross-validations in the repository because the methods obtain nonadiabatic electronic transfer through mathematically different numerical objects.

---

# 20. Incremental PySCF snapshot graph

v0.8 adds a PySCF-facing builder that accepts a newly evaluated CASSCF snapshot and only computes overlap edges to requested existing nodes.

For a new node $v$ and selected neighbors $u$,

$$
\boxed{
O_{uv}
=
\langle\Psi(u)|\Psi(v)\rangle
}
$$

is computed using the full many-electron machinery from v0.6.

The edge is immediately reduced to its polar link and added to the v0.7 graph structure.

This means the graph no longer needs to be rebuilt from scratch after every spawn.

---

# 21. Public raw PySCF point + snapshot API

The v0.6 tracked backend is extended with

```python
evaluate_raw_with_snapshot(geometry)
```

which returns

```text
CartesianElectronicStructurePoint
CASSCFWavefunctionSnapshot
```

without imposing a sequential state label.

That is the correct input for a graph-based calculation because state/gauge consistency is then established by the overlap graph itself.

---

# 22. Cartesian graph mode

A raw PySCF point has

$$
\nabla_RE_I
$$

and

$$
d_{IJ}^{(R)}
$$

in all $3N$ Cartesian coordinates.

The incremental snapshot graph can directly flatten these to

$$
(n_{\rm state},3N)
$$

and

$$
(n_{\rm state},n_{\rm state},3N)
$$

when a Cartesian graph calculation is desired.

Generalized-coordinate projection from v0.5 remains available when a reduced coordinate set is preferred.

---

# 23. Raw-overlap diagnostics remain essential

Before replacing a raw overlap with its polar link, v0.8 records its singular values

$$
\sigma_k(O)
$$

and unitarity defect

$$
\boxed{
\epsilon_U
=\|O^\dagger O-I\|_F.
}
$$

A polar factor is always unitary; that fact must not be allowed to hide a poor underlying electronic-state overlap.

Large singular-value loss means the selected electronic manifold is changing or the geometry step is too large.

---

# 24. What v0.8 demonstrates

v0.8 now has an executable benchmark in which:

- TBF centers move;
- new electronic graph nodes are appended every step;
- centroid nodes are appended for each interacting pair;
- graph cycles accumulate over time;
- a child TBF is dynamically spawned;
- the child enters with zero coefficient;
- coupled propagation gives the child nonzero amplitude;
- $C^\dagger SC$ remains numerically stable.

For the default short benchmark, graph growth is monotonic and the norm remains at unity to essentially floating-point precision.

---

# 25. What v0.8 still does not claim

This remains an intentionally readable research framework.

It is **not yet production AIMS** because it does not include:

1. production optimal-spawning algorithms;
2. a full SPA0/SPA1 saddle-point hierarchy;
3. adaptive graph pruning;
4. large ensemble initial conditions;
5. fully variable Gaussian width matrices in the graph-AIMS runner;
6. distributed ab initio evaluation;
7. asynchronous failure recovery;
8. a global multidimensional graph-gauge optimization updated after every noisy ab initio edge;
9. a rigorous higher-order electronic contribution to $T_{ij}$ beyond the metric-compatible minimal connection.

The scientific advance in v0.8 is narrower:

$$
\boxed{
\text{the gauge graph is now part of the actual time propagation rather than only a static diagnostic.}
}
$$
