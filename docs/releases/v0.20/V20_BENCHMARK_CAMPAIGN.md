# v0.20 Sparse Molecular Benchmark Campaign

## Goal

The campaign tests whether v0.19 molecular data can flow through Gaussian dynamics
without constructing every molecular Gaussian pair and pair centroid.

## 1. Canonical moving sparse-vs-dense propagation

```text
N = 20
possible off-diagonal pairs = 190
dt = 0.002
steps = 20
scrambled raw electronic roots + tracked indexed provider
```

Sparse controls:

```text
enter score = 0.030
exit score = 0.015
search overlap floor = 1e-5
local score budget = 0.010
audit every 5 steps
4 priority + 4 random omitted-edge samples
```

The dense reference uses the exact same molecular pair-centroid approximation but
retains all pairs.

Dense sentinels use separate provider instances.

## 2. Molecular matrix convergence

A fixed 20-Gaussian snapshot is compared with a complete dense molecular reference.

Score-threshold axis:

```text
0.12, 0.08, 0.05, 0.03, 0.02, 0.01, 0.005
```

Local-budget axis at enter score `0.03`:

```text
unbounded, 0.05, 0.02, 0.01, 0.005, 0.0
```

## 3. Bounded-locality scaling

An irregular two-dimensional chain avoids duplicate centroid geometries.

```text
N = 20, 40, 80, 160
spacing ~ 2
small deterministic x/y jitter
```

The campaign records active edges, exact molecular pair checks, new electronic points,
formal dense pairs, KD-tree diagnostics, and sparsity.

## 4. Online controller stress test

A two-Gaussian problem deliberately starts with:

```text
search_overlap_floor = 0.90
```

The sampled audit must detect a hidden important edge, lower the search floor, rebuild,
and pass the immediate re-audit.

## Reproduction

```bash
python examples/102_recompute_v020_campaign.py
```

Canonical output:

```text
results/v020_sparse_molecular_campaign.json
```
