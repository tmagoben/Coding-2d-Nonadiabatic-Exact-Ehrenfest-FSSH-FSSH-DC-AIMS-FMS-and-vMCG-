# v0.21.3 Validation Contract

Machine-readable campaign:

```text
results/v0213_soc_contract_freeze_campaign.json
```

All **20/20** v0.21.3 checks pass, including inherited v0.21.2 acceptance. The
cumulative package suite contains **263 passing tests**.

## Defects closed

The strict fixture has Hermiticity residual

```text
8.944271909998265e-07
```

NumPy's historical default `allclose` relative tolerance accepts that fixture even when
only `atol=1e-12` is supplied. The v0.21.3 explicit residual validator rejects it at the
$10^{-12}$ contract threshold.

At exact two-state degeneracy, the density-guided force is $-1$ in the base frame. A
coordinate-dependent complex gauge gives

```text
current force covariance error:  2.220446049250313e-16
retained force covariance error: 1.1102230246251565e-16
```

After the electronic coefficient block is set exactly to zero, its transported guide
density retains the same physical force. A never-seeded zero block has exactly zero
force. The old lowest-eigenvector fallback is rejected.

## Electronic contract

The campaign declares a four-component model space containing one complete singlet and
one complete triplet. It verifies:

| Gate | Result |
|---|---:|
| zero-SOC H composition error | 0.0 |
| zero-SOC K composition error | 0.0 |
| maximum H/K/D/mass structural residual | 0.0 |
| cm$^{-1}$ to hartree round-trip error | 0.0 |
| incomplete triplet rejected | yes |
| nonzero D in a fixed frame rejected | yes |
| nonzero SOC without SOC provenance rejected | yes |

## Arbitrary-state projection

A one-dimensional, four-state complex target built from a known Gaussian gives

```text
fidelity:                  1.0
relative residual:         1.3638610398823859e-27
coefficient-vector error:  3.6931911946112265e-14
metric condition number:   1.0
```

The regression suite separately exercises a two-dimensional, three-state target.

## Fingerprinted complex cache

For a three-state complex operator fixture:

```text
misses on first provenance: 1
hits on first provenance:   1
misses on changed provenance: 1
base provider calls:        2
complex round-trip error:   0.0
imaginary signal retained:  0.011539731181531185
```

Changing one provenance parameter creates a distinct fingerprint and a distinct cache
entry at the same geometry.

The integrated three-step v0.21.3 runner performs 9 density-guide state rollbacks
between corrector trials, commits 5 accepted coefficient refreshes, and gives maximum
norm drift $5.551115123125783\times10^{-16}$.

## Acceptance gates

The 20 gates cover:

1. strict invariant rejection of the old relative-tolerance defect;
2. complete-multiplet model-space validation;
3. zero-SOC H equivalence;
4. zero-SOC K equivalence;
5. strict operator structure;
6. explicit unit conversion;
7. fixed-frame D semantics;
8. explicit SOC provenance;
9. current-force covariance at degeneracy;
10. retained-force covariance at degeneracy;
11. retained-density force correctness;
12. unseeded zero-block behavior;
13. retirement of the unsafe eigenvector fallback;
14. transactional predictor/corrector guide state;
15. arbitrary-state projection fidelity;
16. arbitrary-state projection residual;
17. exact complex cache round trip;
18. preservation of imaginary operator data;
19. cache separation by provenance;
20. inherited v0.21.2 acceptance.

## Scope boundary

PySCF is not installed in the build environment and no real PySCF trajectory is
claimed. The campaign contains no physical SOC matrix. Passing these gates establishes
interface readiness, not SOC dynamics, ab-initio SOC accuracy, or production AIMS.
