# v0.10 Benchmark Results

This file records the deterministic compact release campaign generated during the
v0.10 build on 2026-08-12.

## Near-CI passage

The release benchmark uses

```text
q0 = (-0.60, 0.25)
p0 = (10.0, 0.0)
A  = diag(1.4, 1.4)
M  = 5
t_final = 0.60
```

The nonzero y impact parameter avoids placing the TBF center exactly on the singular
adiabatic CI point while still producing strong nonadiabatic transfer in the exact
wavepacket.

## Exact versus managed reduced electronic density

The preferred analytic-model observable is the reduced electronic density in the
global diabatic basis.

```text
exact diabatic populations:   [0.22600611046735578, 0.7739938895326441]
managed diabatic populations: [0.9502784027694906, 0.049721597230509486]

exact purity:   0.6762081969769371
managed purity: 0.9748852470781482

exact linear entropy:   0.3237918030230629
managed linear entropy: 0.025114752921851835
```

The managed reference setting therefore remains far from the exact electronic reduced
state for this demanding passage.

This is an intentional and scientifically useful result: v0.10 does not reinterpret
norm conservation as physical accuracy.

## Sensitivity/error budget

```text
total vs exact:              1.0242756986247297
exact discretization proxy:  5.760829006233664e-05
managed timestep proxy:      0.00038484812098345305
SPA truncation proxy:        0.0008869391506525756
spawn-threshold proxy:       0.0
basis-size proxy:            0.5346564941377776
dominant controlled proxy:   basis_size_proxy
```

These quantities are correlated sensitivity probes and are not added in quadrature.

The dominant controlled sensitivity in this compact campaign is the Gaussian basis
size.  The total discrepancy remains substantially larger than the controlled
timestep/SPA/exact-grid proxies, indicating a structural basis/branching deficiency
rather than a simple time-integration error.

## Acceptance result

```json
{
  "checks": {
    "conditioning": true,
    "exact_reference_population": false,
    "norm": true,
    "population_sum": true,
    "pruning_loss": true
  },
  "metrics": {
    "final_norm": 0.9999999194441826,
    "final_populations": [
      0.01298193428853648,
      1.0194141115734356
    ],
    "max_basis_size": 4,
    "max_condition_number": 19927.274562082956,
    "max_norm_error": 8.055581735000317e-08,
    "max_spa1_relative_correction": 2.2487045511537406e-05,
    "observed_population_l2_vs_reference": 1.0242756986247297,
    "observed_populations": [
      0.9502784027694906,
      0.049721597230509486
    ],
    "population_l2_vs_reference": null,
    "prune_count": 0,
    "spawn_count": 3,
    "total_pruning_loss": 0.0
  },
  "passed": false
}
```

The benchmark **fails the exact-reference population criterion** while passing the
configured norm, population normalization, conditioning, and pruning-loss checks.

That separation is exactly the intended v0.10 behavior.

## Exact grid x timestep observation

The exact surface shows extremely small timestep dependence at fixed `grid_n=64`,
while spatial-grid changes are larger for the coordinate-dependent adiabatic
population observable.

Therefore the finest exact row is retained as a **candidate reference**, not claimed
as an analytically exact answer.

The global-diabatic reduced-density population used in the compact error budget is
considerably less sensitive between the selected fine/next-coarse exact calculations,
with the reported population-vector proxy above.

## Machine-readable result

See:

```text
results/v010_compact_release_campaign.json
```
