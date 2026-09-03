# v0.20 Build Validation Report

Validated on 2026-08-13.

## Source validation

```text
330 Python files parsed successfully with Python AST.
```

## Cumulative automated regression suite

```text
231 passed in 27.43 s
```

The suite remains cumulative from v0.1 through v0.20.

New v0.20 coverage includes:

```text
buffered KD-tree nearest-neighbor exactness
indexed molecular-provider tracking integration
all-edge sparse molecular S/H/T equality with dense reference
subquadratic molecular candidate construction on local chains
sparse moving-basis propagation equality with dense propagation when all edges are retained
online sampled-audit geometric search relaxation
complete v0.20 release acceptance
```

## Canonical campaign

Machine-readable output:

```text
results/v020_sparse_molecular_campaign.json
```

Release status:

```text
passed = True
```

Every configured acceptance check passes.

## Canonical sparse molecular propagation

```text
basis size:
20

possible off-diagonal pairs:
190

active edges:
36

average off-diagonal sparsity:
81.05 %

dense-metric coefficient error:
0.000670432070175427

center error:
0.0

momentum error:
0.0

maximum generalized norm drift:
8.881784197001252e-16
```

## Independent dense sentinels

Initial:

```text
S error:
0.002954919263128477

H error:
0.0016775567863378392

T-seed error:
0.01718013917379461
```

Final:

```text
S error:
0.002954918852389199

H error:
0.0016775567668417452

T-seed error:
0.01718013850955588
```

Both sentinels pass.

The sentinel calculations use independent provider/cache instances so they do not
pre-populate the production sparse electronic cache.

## Sampled molecular audits

```text
sampled audit failures:
0

audit checkpoints:
4
```

Every canonical sampled audit passes.

## Electronic-work reduction

```text
sparse production cache misses:
2149

dense reference cache misses:
4239

new electronic-point reduction:
49.30 %
```

Diagnostic wall time:

```text
sparse:
3.945 s

dense:
7.105 s

diagnostic speedup:
1.80 x
```

Wall time is environment-specific; the electronic-point count is the stronger metric.

## Molecular sparse convergence

Score-threshold convergence:

```text
S monotone:
True

H monotone:
True

T monotone:
True
```

Local-score-budget convergence:

```text
S monotone:
True

H monotone:
True

T monotone:
True
```

Zero-budget final snapshot:

```text
S error:
6.153403390687608e-07

H error:
4.528515298636401e-07

T error:
6.984120676640802e-06
```

## Bounded-locality scaling

At `N=160`:

```text
active edges:
152

exact molecular pair checks:
317

formal dense pairs:
12720

pair-check reduction:
97.51 %

matrix sparsity:
98.81 %

new electronic points:
477
```

Fitted exponents:

```text
active edges:
1.117285787300988

exact molecular pair checks:
1.0324491201728092

new electronic points:
1.0213244854837158

formal dense pairs:
2.0213244854837162
```

The bounded-locality molecular pair-decision path is therefore close to linear in this
benchmark while the formal dense pair count is quadratic.

## Indexed electronic cache

Canonical production diagnostics:

```text
nearest queries:
2148

KD-tree queries:
2148

buffer distance checks:
16086

rebuilds:
135

indexed points:
2145

buffered points at end:
4
```

Automated tests verify exact nearest-neighbor agreement against brute force.

## Online geometric controller

The release deliberately tests an over-aggressive search floor.

```text
initial search floor:
0.9

missed-edge score:
0.45067614545809304

initial audit passed:
False

new search floor:
0.45

immediate re-audit passed:
True
```

This validates active error recovery rather than a passive diagnostic only.

## Representative examples executed

```text
examples/95_v020_sparse_canonical.py
examples/96_v020_molecular_sparse_convergence.py
examples/97_v020_sparse_scaling.py
examples/98_v020_indexed_cache.py
examples/99_v020_online_controller.py
examples/100_v020_electronic_work.py
examples/101_v020_pyscf_sparse_status.py
```

`examples/102_recompute_v020_campaign.py` uses the same benchmark function exercised by
the release acceptance test.

## PySCF status

```text
PySCF installed in build environment:
False

real PySCF runtime validated:
False
```

The v0.20 sparse path is connected to the existing raw SA-CASSCF snapshot and
many-electron overlap bridge, but real PySCF runtime validation remains outside this
build environment.

## Scientific limitations

v0.20 finishes the sparse molecular software machinery but does not yet claim:

```text
production AIMS matrix elements
real PySCF sparse trajectory validation
complex-valued molecular derivative-coupling contracts
spin-orbit coupling dynamics
production-scale asynchronous electronic scheduling
```
