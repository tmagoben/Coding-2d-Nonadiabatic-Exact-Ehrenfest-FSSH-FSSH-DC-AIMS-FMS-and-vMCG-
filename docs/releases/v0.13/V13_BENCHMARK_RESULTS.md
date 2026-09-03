# v0.13 Benchmark Results

The machine-readable release campaign is:

```text
results/v013_residual_driven_campaign.json
```

## 1. Candidate dictionary

```text
candidate count: 726
position radius: 1.0
position spacing: 0.2
width scales: [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
momentum offsets: [[0.0, 0.0]]
```

## 2. Residual-selection ladder

| basis size | projection fidelity | relative residual | initial density error | condition number |
|---:|---:|---:|---:|---:|
| 1 | 0.226644940 | 0.773355060 | 0.803042622 | 1 |
| 2 | 0.593851489 | 0.406148511 | 0.051562340 | 19.1725 |
| 3 | 0.720908053 | 0.279091947 | 0.058239963 | 38.3095 |
| 4 | 0.756081846 | 0.243918154 | 0.035653063 | 67.3017 |
| 5 | 0.813153425 | 0.186846575 | 0.046235074 | 76.5341 |
| 6 | 0.837495763 | 0.162504237 | 0.028910339 | 97.7717 |
| 7 | 0.855394582 | 0.144605418 | 0.045534626 | 228.988 |
| 8 | 0.870329432 | 0.129670568 | 0.037019233 | 241.171 |
| 9 | 0.877417170 | 0.122582830 | 0.034136075 | 899.806 |
| 10 | 0.882251454 | 0.117748546 | 0.033619920 | 1172.28 |
| 11 | 0.890252196 | 0.109747804 | 0.032091403 | 3458.02 |

The Hilbert residual decreases monotonically by construction.

The reduced-density error does not decrease monotonically at every intermediate
basis size, which is expected because a nonlinear reduced observable is not identical
to the Hilbert projection norm.

## 3. Pure residual-greedy result

```text
basis size: 11
projection fidelity: 0.8910387734744106
relative residual: 0.10896122652558952
initial density error: 0.03719393002837838
condition number: 4664.580844663738
```

Pure residual greedy gives the smallest wavefunction residual among the two v0.13
11-Gaussian selection modes tested here.

It does not give the smallest initial electronic-density error.

## 4. Density-screened residual reference

The release reference uses the top-30 residual shortlist followed by an initial-density
screen.

```text
basis size: 11
projection fidelity: 0.8902521956060818
relative residual: 0.10974780439391818
initial density error: 0.03209140317550961
initial condition number: 3458.01502834873
```

## 5. Representation-consistent dynamics

```text
projected-state dynamics density error:
0.00011354880287339317

original-target density error:
0.03178630139393256

trace distance:
0.022476309264489125

population L2 error:
0.025521902605714804

Gaussian populations:
[0.2440528208686393, 0.7559471791313606]

exact target populations:
[0.22600611046735578, 0.7739938895326441]

Gaussian purity:
0.663555959251488

exact target purity:
0.6762081969769371

purity error:
0.0126522377254491

coherence:
[0.12753211562743985, -0.0021402345741528663]

exact target coherence:
[0.11414352736213967, -0.001643821331551652]

coherence magnitude error:
0.013394709617860698

coherence phase error / rad:
0.0023799927838891125

Bloch-vector error:
0.0449526185289782

maximum generalized norm drift:
1.0593904686828637e-06

maximum condition number:
3465.8914579773386
```

The projected-state dynamics error remains far smaller than the finite initial
representation error.

## 6. v0.12 comparison

| metric | v0.12 | v0.13 |
|---|---:|---:|
| basis size | 9 | 11 |
| projection fidelity | 0.832276024 | 0.890252196 |
| relative residual | 0.167723976 | 0.109747804 |
| initial density error | 0.035454580 | 0.032091403 |
| projected dynamics error | 0.000290228693 | 0.000113548803 |
| target density error | 0.035000281 | 0.031786301 |
| target population error | 0.028108993 | 0.025521903 |
| coherence phase error / rad | 0.001960749 | 0.002379993 |
| max condition | 2235.29 | 3465.89 |

v0.13 is not a dramatic numerical leap over v0.12 because v0.12 was already close to
the representation limit of a compact bank.

Its key improvement is that the basis is now selected by an explicit residual-reduction
algorithm instead of by a manually specified nine-point layout.

## 7. TDSE-defect enrichment

The reference defect diagnostic gives:

```text
selected candidate:
dq=(0.0, -0.4);dp=(0.0, 0.0);width_scale=4

defect norm before:
0.31502411763651

defect norm after:
0.28652498129723825

relative defect before:
0.030448012071349142

relative defect after:
0.02769348631061861

capture fraction:
0.17274884030759896

predicted squared defect reduction:
0.01714362854505765

actual squared defect reduction:
0.017143629785278974

expanded condition number:
7460.863531968121
```

The candidate is inserted with zero electronic coefficient, so the instantaneous
wavefunction is unchanged.

The predicted and actual squared-defect reductions agree to relative error

```text
7.234298863186708e-08
```

which directly validates the residual-capture algebra.

## 8. Acceptance

```json
{
  "checks": {
    "coherence_phase": true,
    "conditioning": true,
    "defect_gain_prediction": true,
    "defect_reduction": true,
    "initial_density_representation": true,
    "monotone_residual_refinement": true,
    "norm": true,
    "projected_dynamics": true,
    "target_full_density": true,
    "target_population": true
  },
  "defect_prediction_relative_error": 7.234298863186708e-08,
  "passed": true,
  "thresholds": {
    "max_coherence_phase_error": 0.003,
    "max_condition_number": 5000.0,
    "max_defect_prediction_relative_error": 0.005,
    "max_initial_density_error": 0.033,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.0002,
    "max_target_density_error": 0.033,
    "max_target_population_error": 0.03,
    "min_defect_squared_reduction": 1e-08
  }
}
```

All configured v0.13 acceptance criteria pass.

## 9. Scientific interpretation

v0.13 establishes two useful facts.

First, adaptive Gaussian basis construction can be tied directly to a measurable
projection residual instead of relying only on nonadiabatic-coupling thresholds.

Second, the same residual logic extends to the instantaneous Schrödinger defect after
the dynamics begins.

The release stops before implementing a full time-adaptive state machine.  That is
reserved for the next development stage so that trigger cadence, conditioning,
trajectory guidance, and basis removal can be validated separately.
