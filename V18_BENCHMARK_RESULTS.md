# v0.18 Benchmark Results

Canonical machine-readable output:

```text
results/v018_convergence_complete_campaign.json
```

## Canonical state

```text
final basis size:
13

average basis size:
12.525

projected-reference full-wavefunction fidelity:
0.982566093411826

phase-aligned wavefunction L2 error:
0.13232747836123407

nuclear-density L2 error:
0.052341235444456596

nuclear-density total variation:
0.050223497471400924

centroid L2 error:
0.001269391116081437

covariance Frobenius error:
0.01269760151753039

projected reduced-density error:
0.00010573932284646514

original-target reduced-density error:
0.03329249794783041

original-target population error:
0.028073109470748484

original-target coherence phase error:
0.0028907634670896944

norm drift:
1.2434515149761793e-06

maximum condition number:
6509.218903498147
```

## Original-target full wavefunction

Because the 10-Gaussian initial representation is imperfect, direct comparison against
the original exact target includes that initial error.

```text
initial projection fidelity:
0.8822514544600691

Gaussian vs original-target final fidelity:
0.8691787189021798

Gaussian vs original-target final L2:
0.3679740521281518
```

The exact projected-target fidelity remains constant to approximately machine precision,
confirming that this distinction is real rather than bookkeeping.

## Basis ladder

| Nmax | Final N | Fidelity | Full-wavefunction L2 | Density L2 |
|---:|---:|---:|---:|---:|
| 10 | 10 | 0.963780418 | 0.191190020 | 0.108683249 |
| 11 | 11 | 0.968741225 | 0.177501964 | 0.092296771 |
| 12 | 12 | 0.979107229 | 0.144924302 | 0.054475437 |
| 13 | 13 | 0.982566093 | 0.132327478 | 0.052341235 |

Relative L2 improvement, 10 -> 13 Gaussians:

```text
30.79 %
```

## Timestep self-convergence

| Coarse dt | Fine dt | Phase-aligned solution difference |
|---:|---:|---:|
| 0.01 | 0.005 | 0.000610168531085 |
| 0.005 | 0.0025 | 0.000152815765496 |

Observed order:

```text
1.9974143869640382
```

## Sparse-edge-budget axis

| B_local | Full-wavefunction L2 | Graph sparsity | Final H sentinel error |
|---:|---:|---:|---:|
| 0.030 | 0.145733746 | 0.040389 | 0.0034118532 |
| 0.010 | 0.132327478 | 0.013860 | 0.00013521292 |
| 0.000 | 0.132124297 | 0.000000 | 0 |

## Growth-trigger axis

| Enrich threshold | Final N | Enrichment times/steps | Full-wavefunction L2 |
|---:|---:|---|---:|
| 0.050 | 10 | [] | 0.191190020 |
| 0.035 | 12 | [70, 120] | 0.180089293 |
| 0.030 | 13 | [10, 20, 70] | 0.138203917 |
| 0.025 | 13 | [10, 20, 30] | 0.132327478 |
| 0.015 | 13 | [10, 20, 30] | 0.132327478 |

## Audit-cost reduction

```text
v0.17 dense audit pair factorizations:
506

v0.18 dense sentinel pair factorizations:
146

reduction:
71.15 %

normal sampled audits:
6

sampled audit failures:
0
```

## Batched candidate memory

```text
largest unbatched candidate-grid array:
1044800 elements

largest batched candidate-grid array:
25600 elements

peak reduction:
97.55 %
```

## Acceptance

```text
passed = True
```

All configured checks pass.
