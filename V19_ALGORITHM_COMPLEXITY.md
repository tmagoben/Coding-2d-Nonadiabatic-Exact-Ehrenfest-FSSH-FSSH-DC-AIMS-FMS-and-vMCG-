# v0.19 Algorithmic Complexity

v0.19 introduces molecular electronic-structure cost into a Gaussian dynamics code
base. The relevant complexity hierarchy is therefore broader than nuclear Gaussian
pair algebra alone.

## Symbols

| Symbol | Meaning |
|---|---|
| $N$ | Gaussian/TBF count |
| $P$ | active Gaussian pair count |
| $n_q$ | generalized nuclear coordinates |
| $n_R=3N_{atom}$ | Cartesian coordinates |
| $n_s$ | electronic states |
| $m$ | finite electronic representation dimension in validation |
| $N_c$ | cached electronic points |
| $C_{ES}$ | cost of one electronic-structure calculation |

## 1. Coordinate map

Geometry construction

$$
R(q)=R_0+Jq
$$

costs

$$
O(n_Rn_q).
$$

Projecting all state gradients costs approximately

$$
O(n_sn_Rn_q).
$$

Projecting all derivative couplings costs

$$
oxed{
O(n_s^2n_Rn_q).
}
$$

For SA-CASSCF, this linear algebra is negligible relative to electronic-structure work.

## 2. Electronic backend

A cache miss costs

$$
oxed{
C_{ES}
}
$$

plus overlap/state-tracking work.

For the real intended backend this can include

```text
SCF
SA-CASSCF
state gradients
NACs
cross-geometry wavefunction overlap data
```

No universal polynomial in atom count is claimed because practical CASSCF cost depends
strongly on active space, orbital count, integral algorithms, and convergence.

## 3. Exact cache hits

The quantized-coordinate dictionary lookup is approximately

$$
O(1)
$$

average-case.

The release scan makes 17 repeated exact queries and records 17 cache hits.

The scrambled provider ultimately holds

```text
cache size: 23
backend attempts: 23
cache hits: 17
```

after the scan plus center/centroid graph requests.

## 4. Nearest tracking anchor

The current transparent nearest-anchor implementation scans the trusted cache:

$$
oxed{
O(N_cn_q)
}
$$

per new point.

This is intentionally documented as a remaining scaling issue.

For actual ab-initio work, $C_{ES}$ is normally dominant for modest $N_c$, but a
large persistent molecular cache should eventually use a spatial index or structured
trajectory graph rather than an unbounded linear search.

## 5. Cross-geometry state overlap

For the deterministic validation backend with state matrix

$$
V\in\mathbb C^{m\times n_s},
$$

the overlap

$$
O=V_A^\dagger V_B
$$

costs

$$
O(mn_s^2).
$$

For CASSCF snapshots, the many-electron overlap engine includes cross-AO/MO overlap and
CI determinant-space contractions. Its cost can be substantially larger and is not
represented by the finite-vector benchmark.

## 6. State assignment

Legacy small-state tracking enumerated permutations:

$$
O(n_s!).
$$

v0.19 uses the Hungarian algorithm for the best assignment:

$$
O(n_s^3).
$$

To retain the exact second-best assignment margin, it solves $n_s$ constrained
assignments:

$$
oxed{
O(n_s^4).
}
$$

The 16-state diagnostic completes with a valid permutation.

Measured wall time:

```text
0.000302138 s
```

This wall time is diagnostic only.

The complexity change from factorial to polynomial is the portable result.

## 7. Center-centroid electronic graph

A complete basis of $N$ TBFs has

$$
\frac{N(N-1)}2
$$

pair centroids.

The validation graph therefore contains

$$
N+\frac{N(N-1)}2
$$

electronic nodes and

$$
N(N-1)
$$

center-centroid edges.

Thus a complete pair-centroid graph is

$$
oxed{
O(N^2)
}
$$

in both electronic nodes and edges.

That is acceptable for the small v0.19 validation basis, but **not** the intended
large-basis production architecture.

The later molecular sparse implementation should create centroid electronic points only
for active Gaussian locality edges:

$$
P\ll N^2.
$$

Then centroid graph cost becomes approximately

$$
O(N+P).
$$

## 8. Direct local Gaussian matrices

The inherited local pair approximation evaluates all ordered Gaussian pairs:

$$
O(N^2)
$$

pair matrix elements.

Each pair may query a centroid electronic point.

With caching, repeated identical geometry requests avoid repeated backend calls, but
the v0.19 direct runner is still a small-basis reference implementation rather than the
final sparse molecular propagator.

## 9. Electronic scheduling cost model

Candidate estimates are deliberately separated into

```text
cached point  -> normalized cost 0.1
nearby point  -> normalized cost 0.5
new point     -> normalized cost 5.0
```

The current values are dimensionless.

For a PySCF release they should be calibrated from observed timings:

$$
C_{candidate}
=
C_{SCF}
+
C_{CASSCF}
+
C_{grad}
+
C_{NAC}
+
C_{overlap}
$$

as appropriate.

## 10. Failure handling

Retry cost is

$$
O((r+1)C_{ES})
$$

for $r$ configured retries.

The default has no silent fallback.

The optional bounded fallback adds a nearest-cache search

$$
O(N_cn_q)
$$

but no new electronic calculation.

## 11. Complexity summary

The v0.19 architecture has removed one major molecular scaling problem:

$$
oxed{
n_s!
\rightarrow
n_s^4
}
$$

for tracked-state assignment.

The dominant unresolved scaling issues are now:

1. electronic-structure cost $C_{ES}$;
2. linear nearest-anchor cache search for very large caches;
3. complete $O(N^2)$ molecular centroid graphs in the validation builder;
4. inherited $O(N^2)$ local molecular Gaussian matrix construction;
5. many-electron overlap cost for larger CASSCF spaces.

These are the correct issues to address before adding SOC complexity.
