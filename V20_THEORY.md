# v0.20 Theory: End-to-End Sparse Molecular Gaussian Dynamics

v0.20 is the **sparse molecular machinery completion release**.

v0.19 established the molecular electronic-structure and state-tracking bridge.
v0.20 removes the remaining requirement to build every molecular Gaussian pair and
every molecular pair centroid during propagation.

The release architecture is

$$
\boxed{
\text{Gaussian geometry screen}
\rightarrow
\text{molecular S/H/T candidate score}
\rightarrow
\text{local score budget}
\rightarrow
\text{sparse CSR S/H/T}
\rightarrow
\text{sampled omitted-edge audit}
\rightarrow
\text{sparse moving-basis propagation}.
}
$$

The dense molecular pair graph is now a **validation reference**, not the normal
propagation path.

## 1. Sparse molecular TBF graph

Each Gaussian carries a persistent uid,

$$
G_i=(u_i,a_i,q_i,p_i,A_i),
$$

where $a_i$ is the local tracked adiabatic state.

For $N$ basis functions,

$$
P_{\rm dense}=\frac{N(N-1)}2.
$$

v0.20 evaluates molecular pair-centroid electronic structure only for geometrically
local candidates plus fixed-size audit samples.

## 2. Conservative geometric pre-screen

For positive-definite widths define

$$
a_i=\lambda_{\min}(A_i),\qquad
a_j=\lambda_{\min}(A_j),
$$

and

$$
h_{ij}=\frac1{1/a_i+1/a_j}.
$$

The nuclear overlap satisfies

$$
\boxed{
|S_{ij}^{\rm nuc}|
\le
\exp\left[-\frac12 h_{ij}\|q_i-q_j\|^2\right].
}
$$

Momentum mismatch can only reduce the exact Gaussian overlap magnitude.

A KD-tree radius generated from `search_overlap_floor` therefore removes distant
pairs without molecular pair-centroid electronic calculations.

This is a candidate-generation bound on the **nuclear overlap**, not a rigorous bound
on the full molecular $H_{ij}$ or moving-basis $T_{ij}$.

## 3. Pair-centroid electronic transport

For candidate pair $(i,j)$,

$$
q_c=\frac12(q_i+q_j).
$$

The molecular provider supplies accepted tracked electronic snapshots at
$q_i$, $q_j$, and $q_c$.

Let

$$
O_{ci}=\langle\Phi(q_c)|\Phi(q_i)\rangle,
\qquad
O_{cj}=\langle\Phi(q_c)|\Phi(q_j)\rangle.
$$

Their unitary polar factors are

$$
U_{ci}=\operatorname{polar}(O_{ci}),
\qquad
U_{cj}=\operatorname{polar}(O_{cj}).
$$

The local adiabatic unit vectors are transported to the centroid frame:

$$
|v_i^c\rangle=U_{ci}|e_{a_i}\rangle,
\qquad
|v_j^c\rangle=U_{cj}|e_{a_j}\rangle.
$$

Define

$$
s_{ij}^{e}=\langle v_i^c|v_j^c\rangle,
$$

and

$$
v_{ij}^{e}
=
\langle v_i^c|H_e(q_c)|v_j^c\rangle.
$$

These quantities are gauge invariant under consistent local electronic gauge changes.

## 4. Sparse molecular pair matrices

The v0.20 discrete pair-centroid approximation uses

$$
\boxed{
S_{ij}=S_{ij}^{\rm nuc}s_{ij}^{e},
}
$$

$$
\boxed{
H_{ij}
=
T_{ij}^{\rm nuc}s_{ij}^{e}
+
S_{ij}^{\rm nuc}v_{ij}^{e},
}
$$

and the oriented moving-basis seed

$$
\boxed{
T_{ij}^{(0)}
=
\tau_{ij}^{\rm nuc}s_{ij}^{e}.
}
$$

Only diagonal entries and active graph edges are stored in sparse CSR matrices.

This remains the repository's **discrete local-diabatic pair-centroid approximation**.
It is not the full continuous AIMS matrix-element expression.

## 5. Molecular S/H/T edge score

Define

$$
s_{ij}=|S_{ij}|,
$$

$$
h_{ij}
=
\frac{|H_{ij}|}
{\max\left[\sqrt{|H_{ii}H_{jj}|},E_{\rm floor}\right]},
$$

and

$$
t_{ij}
=
\Delta t
\sqrt{
|T_{ij}^{(0)}|^2+
|T_{ji}^{(0)}|^2
}.
$$

The local molecular importance score is

$$
\boxed{
\eta_{ij}
=
\sqrt{
(w_Ss_{ij})^2+
(w_Hh_{ij})^2+
(w_Tt_{ij})^2
}.
}
$$

The canonical weights are

$$
w_S=1,\qquad w_H=0.20,\qquad w_T=1.
$$

New edges require $\eta_{ij}\ge\eta_{\rm enter}$; existing edges persist down to the
smaller exit threshold.

## 6. Global local-omission budget

For locally scored but omitted candidate edges,

$$
\boxed{
B_{\rm local}
=
\left(\sum_{e\in D}\eta_e^2\right)^{1/2}.
}
$$

If this exceeds the configured budget, the largest omitted scores are promoted back
into the active graph until the budget is satisfied.

The canonical v0.20 budget is

$$
B_{\rm local}\le0.01.
$$

This controls accumulation inside the **scored candidate set**. Geometrically screened
pairs are handled by the independent sampled audit layer.

## 7. Sampled omitted-edge audit and controller

Normal propagation does not perform a dense molecular electronic rebuild.

At each audit checkpoint v0.20 scores two omitted-edge samples:

1. priority edges near the geometric search boundary;
2. deterministic pseudo-random omitted edges.

A violation occurs if an omitted sampled edge satisfies

$$
\eta_{ij}>\eta_{\rm enter}.
$$

That means the score would retain the edge, so the failure belongs specifically to the
geometric search layer.

v0.20 responds with

$$
\tau_{\rm search}\leftarrow r\tau_{\rm search},
\qquad 0<r<1,
$$

rebuilds the graph, and immediately re-audits.

The release stress test intentionally starts at

```text
search_overlap_floor = 0.90
```

and finds a missed edge with score

```text
0.450676145
```

before automatically relaxing to

```text
search_overlap_floor = 0.45
```

and passing the re-audit.

## 8. Dense sentinels use independent caches

Dense validation can itself destroy a sparse cost measurement if it fills the
production electronic cache with every pair centroid.

v0.20 therefore supports separate providers for the initial and final dense sentinels.

The canonical architecture is

```text
production sparse provider
+
independent initial dense-sentinel provider
+
independent final dense-sentinel provider
```

so validation does not pre-warm the production cache.

## 9. Sparse moving-basis propagation

For endpoint matrices $S_n,H_n$ and $S_{n+1},H_{n+1}$, average the nuclear/electronic
seed connection,

$$
T_0=\frac12\left(T_n^{(0)}+T_{n+1}^{(0)}\right).
$$

The metric-compatible correction is

$$
\boxed{
T
=
T_0
+
\frac12\left[
\dot S-T_0-T_0^\dagger
\right],
}
$$

with

$$
\dot S\approx\frac{S_{n+1}-S_n}{\Delta t}.
$$

The coefficients use the sparse midpoint/Cayley solve for

$$
iS\dot C=(H-iT)C.
$$

The nuclear centers use active-surface velocity Verlet.

## 10. Canonical dense-versus-sparse propagation

The release benchmark uses

```text
20 Gaussian TBFs
190 possible off-diagonal pairs
36 active molecular edges
dt = 0.002
20 propagation steps
```

Average off-diagonal sparsity is

$$
\boxed{81.05\%}.
$$

Against the dense molecular reference using the same pair-centroid approximation, the
final phase-aligned coefficient error in the dense overlap metric is

$$
\boxed{\epsilon_C=0.00067043207}.
$$

The nuclear center difference is

$$
0.000e+00,
$$

and generalized norm drift is

$$
\boxed{8.882e-16}.
$$

The final dense sentinel gives

$$
\epsilon_S=0.0029549189,
\qquad
\epsilon_H=0.0016775568,
\qquad
\epsilon_T=0.017180139.
$$

## 11. Electronic-work reduction

The production sparse provider records

```text
new electronic points = 2149
```

while the dense reference records

```text
new electronic points = 4239
```

for the same moving nuclear trajectory.

The reduction is

$$
\boxed{49.30\%}.
$$

This includes sampled sparse audits but excludes dense sentinel calculations.

The measured analytic wall-time ratio is approximately `1.80x`; this is
diagnostic only. The electronic-point count is the more portable metric.

## 12. Molecular sparse convergence

### Score threshold

| Enter score | Active edges | S error | H error | T error |
|---:|---:|---:|---:|---:|
| 0.120 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.080 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.050 | 21 | 0.019303621 | 0.0029620109 | 0.052087976 |
| 0.030 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.020 | 23 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.010 | 27 | 0.0070879361 | 0.0020982307 | 0.029257706 |
| 0.005 | 34 | 0.0034262757 | 0.001697729 | 0.017544481 |

All three matrix errors are nonincreasing as the score threshold is lowered.

### Local score budget

| Budget | Active edges | Promoted | Remaining score L2 | S error | H error | T error |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+09 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.05 | 23 | 0 | 0.037174371 | 0.011138901 | 0.0025357232 | 0.042026352 |
| 0.02 | 29 | 6 | 0.019800169 | 0.005908745 | 0.0020486794 | 0.028194466 |
| 0.01 | 36 | 13 | 0.0098800342 | 0.002953481 | 0.0016790767 | 0.017194256 |
| 0.005 | 53 | 30 | 0.0048760583 | 0.0014519635 | 0.00071056245 | 0.010443794 |
| 0 | 85 | 62 | 0 | 6.1534034e-07 | 4.5285153e-07 | 6.9841207e-06 |

At zero local budget every geometrically scored candidate is restored. The remaining
errors are only

$$
\epsilon_S=6.153e-07,
\qquad
\epsilon_H=4.529e-07,
\qquad
\epsilon_T=6.984e-06,
$$

because pairs below the very small geometric search floor remain screened.

## 13. Bounded-locality molecular scaling

The irregular chain benchmark avoids duplicate pair-centroid geometries.

| N | Active edges | Exact molecular pair checks | Dense pairs | Pair-check reduction | New electronic points |
|---:|---:|---:|---:|---:|---:|
| 20 | 15 | 37 | 190 | 80.53% | 57 |
| 40 | 32 | 77 | 780 | 90.13% | 117 |
| 80 | 71 | 157 | 3160 | 95.03% | 237 |
| 160 | 152 | 317 | 12720 | 97.51% | 477 |

Fitted scaling is

$$
E_{\rm active}\sim N^{1.117286},
$$

$$
M_{\rm pair}\sim N^{1.032449},
$$

$$
N_{\rm ES}\sim N^{1.021324},
$$

while the formal dense pair count behaves as

$$
P_{\rm dense}\sim N^{2.021324}.
$$

At $N=160$, only `317` molecular pair scores are evaluated instead of
`12720` dense pairs:

$$
\boxed{97.51\%}
$$

pair-check reduction and

$$
\boxed{98.81\%}
$$

off-diagonal matrix sparsity.

Worst-case strongly overlapping Gaussian ensembles can still become dense.

## 14. Indexed electronic cache

v0.19 searched every trusted cache point to find the nearest tracking anchor.

v0.20 introduces an exact buffered KD-tree.

The immutable tree handles indexed points; a bounded recent buffer is searched
directly. Normal nearest-query work is therefore approximately

$$
O(\log N_c+B n_q)
$$

between rebuilds.

The canonical run records

```text
nearest queries = 2148
KD queries = 2148
rebuilds = 135
buffered points at end = 4
```

Unit tests verify exact nearest-neighbor agreement with brute force.

## 15. What "sparse machinery complete" means

v0.20 now has a coherent sparse molecular chain:

```text
indexed electronic cache
        ↓
Gaussian candidate search
        ↓
molecular S/H/T pair scoring
        ↓
local omission budget
        ↓
sparse S/H/T storage
        ↓
sparse moving-basis propagation
        ↓
sampled omitted-edge audits
        ↓
independent dense validation sentinels
```

The remaining major limitations are now physics and backend-validation issues rather
than an unfinished all-pairs architecture: real PySCF execution, complete AIMS
matrix-element theory, complex electronic/NAC contracts, and SOC.
