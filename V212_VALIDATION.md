# v0.21.2 Validation Contract

Machine-readable campaign:

```text
results/v0212_pre_soc_hardening_campaign.json
```

All **16/16** configured v0.21.2 checks pass, and the inherited v0.21 acceptance campaign
also passes.

## Unequal-width block covariance

```text
S covariance error: 2.311047820503235e-16
H covariance error: 3.1082368759020804e-16
T covariance error: 1.6051159987815854e-16
```

All-edge sparse construction equals the dense unequal-width reference exactly in the
release benchmark:

```text
{'H': 0.0, 'S': 0.0, 'T': 0.0}
```

## Self-consistent complex-gauge dynamics

| dt | steps | gauge-mapped coefficient error | momentum error | max norm drift |
|---:|---:|---:|---:|---:|
| 0.0100 | 5 | 1.946455e-12 | 1.561e-15 | 6.386e-13 |
| 0.0050 | 10 | 4.858956e-13 | 4.565e-16 | 1.592e-13 |
| 0.0025 | 20 | 1.218428e-13 | 1.144e-16 | 4.086e-14 |

Observed orders:

```text
[2.002130949292044, 1.995625423493078]
```

Minimum observed order:

```text
1.995625423493078
```

The nuclear positions agree to the reported numerical precision in every row.

## Generic electronic observable

```text
base expectation: 0.010686571904291245
gauge expectation: 0.010686571904291247
absolute difference: 1.734723475976807e-18
```

The same test deliberately preserves nonzero imaginary content:

```text
max imaginary H: 0.016864365121461722
max imaginary dH: 0.006399344322890202
max imaginary connection: 0.09364391771472179
```

## Subspace-aware provider

```text
subspace checks: 8
ambiguities: 0
minimum singular value: 0.9999999999999991
```

The provider records Procrustes maps but does not force a gauge rotation.

## Adaptive block lifecycle

```text
zero-block birth state change: 0.0
zero-block prune projection loss: 0.0
pruned coefficient error: 0.0
```

## Complex dtype audit

```text
core files scanned: 9
suspicious casts: []
passed: True
```

The three intentional real casts are confined to the inherited spin-free adiabatic
source adapter and are explicitly listed in the campaign JSON.

## PySCF boundary

PySCF is not installed in this build environment and no real v0.21.2 PySCF trajectory is
claimed. This remains an empirical molecular-backend milestone, but it is not required
before adding an **analytic first SOC model** to the now-generalized core.

> Historical note: the v0.21.3 audit subsequently found and closed a gauge-dependent
> low-amplitude fallback and an implicit-relative-tolerance validation gap. Physical SOC
> therefore targets v0.22 through the stricter v0.21.3 contract.
