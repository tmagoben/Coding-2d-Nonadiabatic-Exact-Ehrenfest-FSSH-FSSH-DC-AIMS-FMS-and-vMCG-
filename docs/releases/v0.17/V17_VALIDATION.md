# v0.17 Validation Contract

v0.17 is accepted only if the sparse approximation is controlled at four levels:

1. local edge-score algebra;
2. global local-score budget;
3. dense matrix audits;
4. physical-grid quantum observables.

## 1. Automated unit/regression checks

New tests verify:

```text
safe global overlap radius
pair-specific overlap bound
S/H/T edge scoring
score hysteresis
one-sided score/search relaxation
global local-score L2 promotions
real sparse runner audit/rebuild behavior
threshold convergence
budget convergence
release acceptance logic
```

All older v0.1-v0.16 tests remain cumulative requirements.

## 2. Online audit contract

The release starts intentionally aggressive.

Initial audit:

```text
S relative error = 0.01647878134803449
H relative error = 0.016940324970549453
result = FAIL
```

The controller relaxes once.

Immediate re-audit:

```text
S relative error = 0.0026588687049987378
H relative error = 0.0020460833614795023
result = PASS
```

Thereafter, every scheduled audit through $t=0.6$ passes.

Required final limits:

```text
S relative error     <= 0.006
H relative error     <= 0.006
Snuc relative error  <= 0.006
unresolved audits    = 0
```

## 3. Local importance budget

The configured local score proxy is

```text
B_local <= 0.08
```

Maximum recorded trajectory value:

```text
0.024083537794312975
```

This check does not replace the dense matrix audit.

## 4. Physical accuracy

Release thresholds:

```text
initial density representation      <= 0.035
projected-state dynamics error      <= 0.001
target density error                <= 0.035
target population error             <= 0.03
coherence phase error               <= 0.0035 rad
norm drift                           <= 1e-4
condition number                     <= 5e3
```

Measured:

```text
initial density error:
0.033619920355630904

projected dynamics error:
0.00013361460054442858

target density error:
0.03333954068459557

population error:
0.02819941365898425

coherence phase error:
0.0029095064228609707

norm drift:
2.0053154399235495e-06

maximum condition:
1431.0606683729504
```

## 5. v0.16 regression

The final v0.17 reduced-density matrix differs from the v0.16 release result by

```text
7.09620610556202e-15
```

with a release limit of `5e-4`.

The measured difference is essentially floating-point noise because the online
controller relaxes to the previously validated sparse representation.

## 6. Threshold convergence

The final snapshot is audited at descending edge-score thresholds.

Required:

```text
S error nonincreasing
H error nonincreasing
finest S error <= 0.001
finest H error <= 0.001
```

Measured finest errors:

```text
S: 0.00026621461118036714
H: 0.00015704048353805722
```

## 7. Local-budget convergence

At fixed nominal edge score, the omitted-score budget is tightened.

Required:

```text
S error nonincreasing
H error nonincreasing
```

At zero budget, all locally scored edges are restored and the final snapshot gives zero
dense/sparse S/H difference.

## 8. Scaling contract

At $N=160$ on the bounded-locality chain:

```text
pair-factorization reduction >= 90%
```

Measured:

```text
93.866 %
```

Fitted limits:

```text
active-edge exponent       <= 1.15
exact-score-check exponent <= 1.20
dense-pair exponent        >= 1.90
```

Measured:

```text
active edge: 1.0324491201728092
exact checks: 1.0557080719105292
dense pairs: 1.97980963642741
```

## 9. What is not claimed

The validation does not establish a universal theorem that the same numerical edge
score or audit tolerance is appropriate for every Hamiltonian.

A new molecular/SOC benchmark must recalibrate the score weights and matrix-error
budget against appropriate dense or independently converged references.
