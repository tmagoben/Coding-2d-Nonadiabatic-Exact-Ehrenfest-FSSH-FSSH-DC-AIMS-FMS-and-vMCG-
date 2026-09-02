# v0.14 Validation Contract

v0.14 is accepted only if the new adaptive controller behaves as an audited numerical
error-control mechanism rather than as an uncontrolled basis-growth heuristic.

## 1. Cumulative regression

Every retained v0.1-v0.13 test must remain passing.

The v0.14 suite additionally checks:

- Hermitian half-build equivalence to the full ordered-pair matrix builder;
- exact pair-count reduction;
- asymptotic-complexity metadata;
- exact leave-one-out pruning loss;
- vectorized dynamic defect ranking versus the slower reference score;
- energy conservation of dynamically generated candidate TBF centers;
- control hysteresis validation;
- actual adaptive enrichment in a short propagation;
- norm conservation after adaptive basis insertion;
- complexity-ledger counters;
- release acceptance logic for adaptation, complexity, and pruning.

## 2. Hermitian half-build equality

For the same mixed-width, mixed-guidance-state Gaussian basis,

```text
build_spinor_complete_lvc_matrices(...)
```

and

```text
build_spinor_complete_lvc_matrices_symmetric(...)
```

must produce numerically identical:

```text
S
H
Snuc
```

to tight floating-point tolerance.

The optimization is rejected if it changes the physics.

## 3. Pair-evaluation accounting

For $N$ Gaussian functions:

```text
ordered builder: N^2
v0.14 half builder: N(N+1)/2
```

The release stores both the actual half-build pair count and the ordered-pair
equivalent over the full variable-basis trajectory.

## 4. Dynamic candidate generation

The candidate generator reuses the local energy-conserving placement equations.

Every generated candidate used by the release must satisfy the configured classical
energy tolerance.

Residual ranking is applied only after this physical candidate filter.

## 5. Dynamic candidate ranking

The vectorized residual ranking is checked against the direct v0.13
`defect_candidate_capture(...)` reference on a controlled candidate set.

The best candidate's capture fraction must agree within the configured numerical
tolerance.

This guards against errors in the matrix orientation/conjugation of the vectorized
$K\times G$ contractions.

## 6. Adaptive enrichment

A release enrichment must satisfy all of the following:

```text
relative TDSE defect >= enrichment threshold
basis size < hard limit after any replacement pruning
candidate capture fraction >= minimum
candidate expanded condition <= limit
zero electronic coefficient at insertion
defect after insertion < defect before insertion
```

The final condition is measured by recomputing the TDSE defect in the enlarged
Galerkin space.

## 7. Hysteresis

The settings object requires:

$$
\eta_{\rm add}>\eta_{\rm remove}.
$$

Configurations that reverse or eliminate this ordering are rejected.

## 8. Cooldown

Ordinary growth/pruning decisions must be separated by the configured minimum number of
time steps.

This is an algorithmic stability constraint, not a physical approximation.

## 9. Exact pruning score

For every Gaussian $j$ the leave-one-out loss

$$
L_j
=
\frac{\sum_a|C_{ja}|^2}{(S^{-1})_{jj}}
$$

is tested against an independent direct projection into the basis with Gaussian $j$
removed.

The two values must agree to numerical precision.

## 10. Pruning stress test

The release inserts a nearly redundant Gaussian with zero coefficient.

Pruning must:

```text
remove that unprotected Gaussian
produce negligible represented-wavefunction loss
reduce the overlap condition number
```

This provides a deterministic test of both the coefficient and nonorthogonal geometry
parts of the pruning rule.

## 11. Representation-consistent accuracy

The exact TDSE and adaptive Gaussian calculation begin from the same 10-Gaussian
projected initial state when computing the projected-state dynamics error.

The release separately reports:

```text
initial representation error
projected-state dynamics error
original-target density error
```

as established in v0.12.

## 12. Release thresholds

```text
initial reduced-density error          <= 0.035
projected-state dynamics error         <= 0.003
original-target density error          <= 0.035
population L2 error                    <= 0.03
coherence phase error                  <= 0.0035 rad
generalized norm drift                 <= 1e-4
maximum condition number               <= 5e3
defect enrichments                     >= 1
every enrichment lowers defect         PASS required
Hermitian pair-evaluation reduction    >= 40%
pruning stress projection loss         <= 1e-10
pruning stress improves conditioning   PASS required
```

These are regression thresholds for the analytic benchmark.

They are not universal chemical-accuracy criteria.

## 13. Complexity validation

Complexity reporting is part of the release output.

The ledger must include:

```text
matrix_build_calls
pair_matrix_evaluations
ordered_pair_equivalent
time_matrix_calls
cayley_solve_calls
defect_evaluations
candidate_ranking_calls
candidate_count_scored
enrichment_events
pruning_audits
pruning_events
peak_basis_size
peak_electronic_dimension
peak_candidate_count
timing categories
symbolic asymptotic scaling
```

See `V14_ALGORITHM_COMPLEXITY.md`.

## 14. PySCF scope

The time-adaptive controller is validated on the analytic two-dimensional LVC model.

The inherited PySCF electronic-structure and gauge-tracking infrastructure remains in
the package, but a full molecular TDSE defect is not claimed.

See `V14_PYSCF_ADAPTIVE_BRIDGE.md`.
