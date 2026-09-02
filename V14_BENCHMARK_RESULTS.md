# v0.14 Benchmark Results

The canonical machine-readable campaign is:

```text
results/v014_time_adaptive_defect_campaign.json
```

## 1. Initial basis

The run starts with a residual-selected 10-Gaussian bank:

```text
projection fidelity: 0.8822514544600691
relative wavefunction residual: 0.11774854553993085
initial reduced-density error: 0.033619920355630904
```

Unlike v0.13, the eleventh Gaussian is not present from $t=0$.

## 2. Adaptive event

At step 10, time 0.05, the TDSE defect triggered enrichment.

```text
basis before: 10
basis after: 11

parent uid: 9
new uid: 10
guidance state: 0

candidate:
parent=9;target=0;pos=target_force:0.06;mom=nac;width=1.35

candidate count searched:
560

relative defect before:
0.03238019095782191

relative defect after:
0.029987957971150208

predicted capture fraction:
0.15067690241980206

orthogonal candidate norm:
0.021818634821044958

expanded condition number:
1454.4566005006334
```

The new two-component coefficient is exactly zero at insertion.

## 3. Final accuracy

```text
initial basis size: 10
final basis size: 11
time-average basis size: 10.925

projected-state dynamics density error:
9.527804623132635e-05

original-target density error:
0.03330494031479218

trace distance:
0.02355014914360283

population L2 error:
0.028084897912098693

Gaussian populations:
[0.24586513222993275, 0.7541348677700673]

exact target populations:
[0.22600611046735578, 0.7739938895326441]

Gaussian purity:
0.6613299727370953

exact target purity:
0.6762081969769371

purity error:
0.014878224239841864

coherence phase error / rad:
0.0028906431794148953

maximum generalized norm drift:
2.115581487549534e-06

maximum condition number:
1470.7558920505405
```

The projected-state dynamics error remains much smaller than the finite initial
representation error.

## 4. v0.13 comparison

| metric | v0.13 | v0.14 |
|---|---:|---:|
| initial/static basis size | 11 | 10 |
| final basis size | 11 | 11 |
| average basis size | 11 | 10.925 |
| projection fidelity | 0.890252196 | 0.882251454 |
| initial density error | 0.032091403 | 0.033619920 |
| projected dynamics error | 0.000113548803 | 9.52780462e-05 |
| target density error | 0.031786301 | 0.033304940 |
| population error | 0.025521903 | 0.028084898 |
| coherence phase error | 0.002379993 | 0.002890643 |

v0.14 trades a very small amount of static-basis accuracy for an actual time-dependent
error-control decision.

That tradeoff is explicit rather than hidden.

## 5. Pruning stress test

```text
removed uid:
999999

fractional projection loss:
0.0

condition before:
379319.12346481933

condition after:
67.30166373596352

condition improvement factor:
5636.103216600353
```

The redundant zero-amplitude Gaussian was removed with zero represented-state loss.

## 6. Complexity ledger

```text
total adaptive runtime / s:
11.28900403300031

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

defect evaluations:
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

The detailed time/memory scaling audit is in `V14_ALGORITHM_COMPLEXITY.md`.

## 7. Acceptance

```json
{
  "checks": {
    "adaptive_enrichment": true,
    "coherence_phase": true,
    "conditioning": true,
    "enrichment_reduces_defect": true,
    "hermitian_pair_reduction": true,
    "initial_density_representation": true,
    "low_loss_pruning": true,
    "norm": true,
    "projected_dynamics": true,
    "pruning_improves_condition": true,
    "target_density": true,
    "target_population": true
  },
  "pair_evaluation_reduction": 0.45420136260408783,
  "passed": true,
  "thresholds": {
    "max_coherence_phase_error": 0.0035,
    "max_condition_number": 5000.0,
    "max_initial_density_error": 0.035,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.003,
    "max_pruning_stress_loss": 1e-10,
    "max_target_density_error": 0.035,
    "max_target_population_error": 0.03,
    "min_enrichment_events": 1,
    "min_pair_evaluation_reduction": 0.4
  }
}
```

All configured v0.14 release criteria pass.
