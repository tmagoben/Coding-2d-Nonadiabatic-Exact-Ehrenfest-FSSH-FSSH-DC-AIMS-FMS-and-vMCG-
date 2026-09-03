# v0.14 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 210 Python files in the repository parsed successfully with Python's AST parser.

## Automated regression suite

```text
167 passed in 8.46 s
```

The cumulative suite retains all v0.1-v0.13 tests and adds v0.14 checks for:

- Hermitian half-build equality with the full ordered-pair S/H builder;
- exact pair-evaluation counting;
- symbolic complexity metadata and scaling proxies;
- exact leave-one-out pruning loss versus direct basis projection;
- energy-conserving dynamic defect candidates;
- vectorized defect-candidate capture versus the direct v0.13 reference;
- hysteretic controller validation;
- short-run real defect-triggered enrichment;
- generalized norm conservation across adaptive growth;
- complexity-ledger counters;
- release acceptance checks for defect reduction, pair-work reduction, and pruning.

## Canonical release campaign

The complete v0.14 campaign was recomputed and saved as:

```text
results/v014_time_adaptive_defect_campaign.json
```

All configured acceptance checks pass.

## Adaptive basis event

```text
step: 10
time: 0.05

basis:
10 -> 11

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

candidate count:
560

relative TDSE defect:
0.03238019095782191
->
0.029987957971150208

predicted capture fraction:
0.15067690241980206

zero-coefficient insertion:
True
```

The measured defect decreases after the basis enlargement.

## Final representation-consistent result

```text
initial basis size: 10
final basis size: 11
average basis size: 10.925

projection fidelity:
0.8822514544600691

initial reduced-density error:
0.033619920355630904

projected-state dynamics density error:
9.527804623132635e-05

original-target density error:
0.03330494031479218

population L2 error:
0.028084897912098693

coherence phase error / rad:
0.0028906431794148953

generalized norm drift:
2.115581487549534e-06

maximum condition number:
1470.7558920505405
```

The projected-state dynamics error remains below the original-target error by more
than two orders of magnitude.

## Pruning stress test

```text
removed uid:
999999

fractional represented-wavefunction loss:
0.0

condition number:
379319.12346481933
->
67.30166373596352

condition improvement factor:
5636.103216600353
```

The redundant zero-amplitude Gaussian is removed with exactly zero measured projection
loss.

## Complexity audit

```text
matrix build calls:
122

Hermitian pair evaluations:
7931

ordered-pair equivalent:
14531

pair-evaluation reduction:
45.420 %

moving-basis T builds:
120

Cayley solves:
120

TDSE-defect evaluations:
13

candidate searches:
1

candidates scored:
560

peak basis size:
11

peak electronic dimension:
22
```

Observed campaign timings:

```text
total adaptive run:
11.289004 s

S/H matrix builds:
4.885241 s

moving-basis T matrices:
4.421099 s

TDSE-defect checks:
0.680050 s

candidate ranking:
0.161364 s

Cayley solves:
0.007415 s

pruning audits:
0.002369 s
```

Wall-time values are environment dependent.  The symbolic complexity and exact
operation/pair counts are the portable complexity statements.

The present small-N run is dominated by repeated unequal-width Gaussian pair algebra.
The dense Cayley solve is asymptotically cubic but numerically tiny at electronic
dimension 22.

See:

```text
V14_ALGORITHM_COMPLEXITY.md
```

for the complete time/memory scaling derivation.

## Representative examples

The following examples were executed successfully during final validation:

```text
examples/53_v014_adaptive_event.py
examples/54_v014_complexity.py
examples/55_v014_pruning_stress.py
examples/56_v013_v014_comparison.py
examples/57_v014_pair_count.py
```

The canonical campaign itself was recomputed directly during the build.
`examples/58_recompute_v014_campaign.py` is the user-facing wrapper for that same
workflow.

## Scientific scope

v0.14 is a time-adaptive analytic-LVC Gaussian reference implementation.

The inherited PySCF SA-CASSCF, many-electron overlap, state-tracking, and gauge-graph
layers remain in the cumulative repository, but the release does not claim a
full-dimensional molecular TDSE-defect controller.

See `V14_PYSCF_ADAPTIVE_BRIDGE.md`.
