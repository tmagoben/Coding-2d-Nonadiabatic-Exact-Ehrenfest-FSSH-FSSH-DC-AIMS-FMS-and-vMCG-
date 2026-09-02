# v0.19 Theory: Molecular Direct-Dynamics Integration

v0.19 is the **molecular/direct-dynamics integration release**.

The goal is not to add spin-orbit coupling yet.

The goal is to establish a defensible bridge

$$
\boxed{
\text{molecular geometry}
\rightarrow
\text{electronic structure}
\rightarrow
\text{state/gauge tracking}
\rightarrow
\text{generalized nuclear coordinates}
\rightarrow
\text{Gaussian dynamics}
}
$$

without making electronic state identity depend on arbitrary backend root ordering or
the order in which branched Gaussian centers happen to request calculations.

The release adds:

1. overlap-capable molecular electronic snapshots;
2. nearest-anchor state/gauge tracking for branched center/centroid queries;
3. deterministic geometry caching, cost estimates, and explicit failure policy;
4. scalable maximum-overlap root assignment;
5. a center-centroid electronic gauge graph for molecular Gaussian pairs;
6. a provider-neutral direct-dynamics wrapper;
7. a raw PySCF SA-CASSCF snapshot bridge;
8. a deterministic Cartesian molecular validation backend.

---

## 1. Cartesian molecular data contract

A molecular backend returns

$$
E_I(R),
$$

Cartesian gradients

$$
\nabla_R E_I,
$$

and derivative couplings

$$
d_{IJ}^{(R)}
=
\langle\Phi_I(R)|\nabla_R\Phi_J(R)\rangle.
$$

For a linear generalized-coordinate map

$$
R(q)
=
R_0+Jq,
$$

the projected quantities are

$$
\boxed{
\nabla_q E_I
=
J^T\nabla_R E_I
}
$$

and

$$
\boxed{
d_{IJ}^{(q)}
=
J^T d_{IJ}^{(R)}.
}
$$

The generalized mass matrix is

$$
\boxed{
M_q
=
J^T M_R J.
}
$$

The v0.19 validation model embeds the exact 2D LVC problem into a two-atom Cartesian
geometry with two orthonormal collective modes.

This backend is **molecular-style**, not ab initio.

---

## 2. Why raw electronic root order is insufficient

An electronic-structure code can change root order or electronic phase between nearby
geometries.

If raw roots are used directly, then a nominal state label $I$ need not represent the
same electronic state at neighboring points.

The v0.19 scrambled backend deliberately introduces:

```text
root permutations
real sign flips
```

as functions of geometry.

Without tracking, the maximum state-energy ordering error in the 17-point scan is

$$
\boxed{
0.0354263241856.
}
$$

After overlap tracking, the maximum errors are

$$
\boxed{
\epsilon_E
=
5.204e-18,
}
$$

$$
\boxed{
\epsilon_{\nabla E}
=
5.204e-18,
}
$$

and

$$
\boxed{
\epsilon_d
=
4.441e-16.
}
$$

Thus the deliberately corrupted raw root convention is removed to numerical precision.

---

## 3. Cross-geometry overlap tracking

Let

$$
O_{ij}
=
\langle
\Psi_i^{(ref)}
|
\Psi_j^{(new)}
\rangle.
$$

The tracked assignment maximizes

$$
\boxed{
\sum_i
|O_{i,\pi(i)}|^2.
}
$$

For a selected state mapping $\pi$, the real-gauge phase correction is chosen so the
assigned overlap is positive.

For a general complex gauge,

$$
p_i
=
\frac{
O_{i,\pi(i)}^*
}{
|O_{i,\pi(i)}|
}.
$$

State properties transform as

$$
E_i'
=
E_{\pi(i)},
$$

$$
G_i'
=
G_{\pi(i)},
$$

and

$$
\boxed{
d_{ij}'
=
p_i^*p_j
d_{\pi(i)\pi(j)}.
}
$$

v0.19 remains a **real-state dynamics contract**. Complex tracking phases are
architecturally available in the assignment layer, but a complex electronic/NAC
contract is intentionally deferred until the SOC release.

---

## 4. Branched queries require a different tracking rule

The older v0.6 tracked PySCF backend is sequential:

```text
geometry n
    tracked against
geometry n-1
```

That is appropriate for one geometry scan or one trajectory center.

It is not sufficient for a branched Gaussian basis, because electronic requests may
come from

```text
TBF center 1
TBF center 2
pair centroid 1-2
TBF center 3
pair centroid 1-3
...
```

and call order is not physical time order.

v0.19 therefore aligns every new point to the **nearest accepted cached electronic
anchor** in generalized-coordinate space.

The first accepted geometry defines the root labels and gauge.

After that seed, the benchmark evaluates the scan in a deliberately shuffled order and
still obtains

$$
\epsilon_E
=
5.204e-18,
$$

$$
\epsilon_d
=
4.441e-16.
$$

This is **order-tolerant local tracking**, not a theorem of global path independence.

Around nontrivial Berry/Wilson holonomy, a globally trivial electronic gauge is still
impossible.

---

## 5. Scalable state assignment

Earlier releases used explicit permutation enumeration.

For $n_s$ states this scales as

$$
O(n_s!).
$$

That is acceptable for two or three states but unsuitable as a molecular default.

v0.19 formulates root assignment as a maximum-weight bipartite matching with weights

$$
w_{ij}
=
|O_{ij}|^2.
$$

The best assignment is obtained with the Hungarian algorithm:

$$
O(n_s^3).
$$

To preserve the previous ambiguity diagnostic, v0.19 also computes the exact
second-best global assignment.

Every distinct assignment must omit at least one edge selected by the best assignment.

Therefore each selected best edge is forbidden in turn and the optimal constrained
matching is recomputed.

This gives total complexity approximately

$$
\boxed{
O(n_s^4)
}
$$

while preserving an exact second-best score margin.

The release validates a 16-state assignment:

```text
nstate:
16

valid permutation:
True

diagnostic wall time:
0.000302138 s
```

Wall time is environment dependent; the algorithmic reduction from factorial to
polynomial complexity is the important result.

---

## 6. Center-centroid molecular gauge graph

For Gaussian centers $i$ and $j$, v0.19 constructs an electronic node at the pair
centroid

$$
q_{ij}
=
\frac12(q_i+q_j).
$$

The graph contains links

$$
i
\leftrightarrow
q_{ij}
\leftrightarrow
j.
$$

Electronic states from each TBF center are transported into the common centroid frame.

The discrete local-diabatic Gaussian pair approximation then uses

$$
S_{ij}
\approx
S_{ij}^{nuc}
\langle e_i|e_j\rangle_{q_{ij}},
$$

and

$$
H_{ij}
\approx
T_{ij}^{nuc}
\langle e_i|e_j\rangle
+
S_{ij}^{nuc}
\langle e_i|H_e(q_{ij})|e_j\rangle.
$$

For the 3-Gaussian molecular validation basis, clean and deliberately scrambled raw
electronic calculations give

$$
\boxed{
\|S_{clean}-S_{scrambled}\|_F
=
0.000e+00
}
$$

and

$$
\boxed{
\|H_{clean}-H_{scrambled}\|_F
=
0.000e+00.
}
$$

The Hamiltonian Hermiticity error is only

$$
1.267e-18.
$$

---

## 7. Molecular direct dynamics

The release direct-dynamics wrapper uses the inherited provider-neutral local Gaussian
approximation.

The deterministic benchmark begins on the upper adiabatic surface and produces one
energy-conserving child at

```text
step = 1
time = 0.02
target state = 0
```

Three calculations are compared:

```text
direct analytic generalized LVC provider
clean Cartesian molecular provider
scrambled-root Cartesian molecular provider + v0.19 tracking
```

The final coefficient difference between analytic and scrambled molecular dynamics is

$$
\boxed{
3.301e-16.
}
$$

The TBF-center difference is

$$
4.475e-16,
$$

and generalized norm drift is

$$
\boxed{
7.772e-16.
}
$$

This validates the molecular provider bridge without claiming that the inherited local
matrix-element approximation is full AIMS.

---

## 8. Electronic cache and cost model

An exact repeated geometry returns a cache hit.

A candidate can be classified as

```text
exact cached point
near an existing electronic point
new electronic point
```

with deterministic normalized costs.

The release demonstration is

```text
cached: 0.1
nearby: 0.5
new:    5.0
```

These are **dimensionless scheduling costs**, not measured SA-CASSCF wall times.

A real molecular workflow should calibrate them against observed SCF/CASSCF/gradient/NAC
costs.

---

## 9. Failure policy

Backend failure is explicit.

The default policy is

```text
raise
```

and does not silently continue.

v0.19 also implements an opt-in

```text
nearest_cache
```

fallback subject to a strict maximum geometry distance.

The validation failure occurs at a distance

$$
0.040,
$$

and exactly one fallback is recorded.

A fallback point is **not promoted into the electronic gauge cache as a new trusted
anchor**.

This prevents one stale approximation from silently contaminating future state
tracking.

---

## 10. PySCF bridge

v0.19 introduces

```text
PySCFRawSnapshotBackendV19
```

which exposes the existing raw SA-CASSCF calculation together with the many-electron
CASSCF wavefunction snapshot.

Cross-geometry overlaps use

```text
casscf_state_overlap_matrix
```

from the v0.6 machinery.

The internal derivative-coupling convention remains

$$
d_{ij}
=
\langle\Phi_i|\nabla_R\Phi_j\rangle.
$$

PySCF's state tuple is still interpreted as

```text
state = (ket, bra)
```

so internal $d_{ij}$ requests use

```text
state = (j, i)
mult_ediff = False
```

for dynamics.

PySCF is **not installed in the v0.19 build environment**.

Therefore the bridge API and overlap-engine path are regression tested, but no real
PySCF runtime result is claimed.

---

## 11. Scope

v0.19 establishes the molecular software and gauge architecture needed before SOC.

It does **not** yet claim:

- production molecular AIMS;
- full continuous AIMS matrix elements;
- globally path-independent state tracking;
- asynchronous electronic-structure scheduling;
- real PySCF runtime validation;
- complex spinor/SOC dynamics.

Those distinctions are deliberate.
