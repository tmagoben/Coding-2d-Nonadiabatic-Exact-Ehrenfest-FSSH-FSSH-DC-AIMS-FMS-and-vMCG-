# v0.17 Benchmark Results

Canonical output:

```text
results/v017_sparse_error_control_campaign.json
```

## 1. Physical result

```text
initial basis:
10

final basis:
11

average basis:
10.925

projection fidelity:
0.8822514544600691

initial density error:
0.033619920355630904

projected-state dynamics error:
0.00013361460054442858

target density error:
0.03333954068459557

target population error:
0.02819941365898425

coherence phase error:
0.0029095064228609707

norm drift:
2.0053154399235495e-06

maximum condition:
1431.0606683729504
```

## 2. Online audit controller

Initial audit:

```text
S error:
0.01647878134803449

H error:
0.016940324970549453

status:
FAIL
```

Automatic relaxation:

```text
enter score:
0.060 -> 0.030

exit score:
0.030 -> 0.015

search floor:
1e-5 -> 5e-6
```

Immediate re-audit:

```text
S error:
0.0026588687049987378

H error:
0.0020460833614795023

status:
PASS
```

Every later audit at steps 20, 40, 60, 80, 100, and 120 passes.

## 3. Final audit

```text
S relative error:
0.005191742661052565

H relative error:
0.003962632349871911

Snuc relative error:
0.005191742661052565

omitted off-diagonal pairs:
3

maximum omitted overlap:
0.02057350476995086

maximum omitted H block:
0.2318916964441307
```

## 4. Local score budget

Configured:

```text
B_local <= 0.08
```

Maximum recorded trajectory value:

```text
0.024083537794312975
```

No trajectory record violates the local proxy budget.

## 5. Physical comparison with v0.16

```text
final reduced-density difference:
7.09620610556202e-15
```

The result is equal to the accepted v0.16 trajectory to floating-point precision.

## 6. TDSE-defect basis event

```text
step:
10

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

basis:
10 -> 11

relative defect:
0.03255613084426919
->
0.03027005382147624

capture fraction:
0.14896413884582896

zero coefficient insertion:
True
```

## 7. Threshold convergence

| Enter score | Active edges | Omitted score L2 | S error | H error |
|---:|---:|---:|---:|---:|
| 0.120 | 46 | 0.19361725 | 0.042094541 | 0.039432122 |
| 0.080 | 48 | 0.11623692 | 0.025222344 | 0.024676408 |
| 0.060 | 49 | 0.093505195 | 0.02029125 | 0.019847861 |
| 0.040 | 51 | 0.042269048 | 0.0092323874 | 0.0073814687 |
| 0.030 | 52 | 0.023742025 | 0.0051917427 | 0.0039626323 |
| 0.020 | 53 | 0.011080847 | 0.0024218886 | 0.001855452 |
| 0.010 | 54 | 0.001212326 | 0.00026621461 | 0.00015704048 |

The S/H errors decrease monotonically as the graph is relaxed.

## 8. Local-budget convergence

| Budget | Active edges | Promoted | Omitted score L2 | S error | H error |
|---:|---:|---:|---:|---:|---:|
| 1e+30 | 49 | 0 | 0.093505195 | 0.02029125 | 0.019847861 |
| 0.1 | 49 | 0 | 0.093505195 | 0.02029125 | 0.019847861 |
| 0.08 | 50 | 1 | 0.072297083 | 0.015669953 | 0.015701252 |
| 0.05 | 51 | 2 | 0.042269048 | 0.0092323874 | 0.0073814687 |
| 0.03 | 52 | 3 | 0.023742025 | 0.0051917427 | 0.0039626323 |
| 0.01 | 54 | 5 | 0.001212326 | 0.00026621461 | 0.00015704048 |
| 0 | 55 | 6 | 0 | 0 | 0 |

At zero local budget all locally scored pair blocks are restored and the final snapshot
matches the dense S/H matrices exactly.

## 9. Bounded-locality construction scaling

At `N=160`:

```text
active edges:
317

all off-diagonal pairs:
12720

exact S/H/T score checks:
630

pair factorizations:
790

dense canonical pairs:
12880

pair reduction:
93.866 %

sparse H density:
0.03099609375

dense assembly:
3.823627 s

sparse assembly:
0.029290 s

diagnostic assembly speedup:
130.54 x
```

Fitted exponents:

```text
active edges:
1.0324491201728092

spatial candidates:
1.0557080719105292

exact S/H/T checks:
1.0557080719105292

pair factorizations:
1.043904194562303

dense pair count:
1.97980963642741
```

## 10. Audit overhead

```text
dense audits:
8

audit pair factorizations:
506

audit time:
0.148051 s
```

This overhead is deliberate and separately accounted.

## 11. Acceptance

All configured v0.17 checks pass.
