# v0.12 Benchmark Results

The machine-readable release campaign is

```text
results/v012_representation_consistent_campaign.json
```

## 1. Representation audit at time zero

The exact strong-CI benchmark starts from a coordinate-dependent adiabatic packet.

Its initial reduced electronic purity is

```text
0.6764597760317345
```

A one-spinor center-frozen replacement has

```text
purity: 0.9999999999999998
initial reduced-density error: 0.28703562527170995
center-frozen populations: [0.03846153846153845, 0.9615384615384617]
```

Thus a substantial part of the older full-density discrepancy exists before
propagation.

## 2. Projection ladder

| bank | Gaussians | projection fidelity | initial density error | projected-dynamics error | final target density error | target population error | condition |
|---|---:|---:|---:|---:|---:|---:|---:|
| one | 1 | 0.226644940 | 0.803042622 | 7.33498788e-05 | 0.802920863 | 0.624351962 | 1 |
| five | 5 | 0.803744614 | 0.0656984625 | 0.00042199215 | 0.0649948568 | 0.0577538825 | 779.33 |
| nine | 9 | 0.832276024 | 0.0354545799 | 0.000290228693 | 0.0350002807 | 0.028108993 | 2235.29 |

The key trend is that the projected-state propagation error is already tiny for every
bank, while the error relative to the original target follows the quality of the
initial representation.

## 3. Nine-Gaussian reference

```text
initial wavefunction fidelity:             0.832276023595292
initial relative wavefunction residual:    0.16772397640470793
initial reduced-density error:              0.03545457994295867

exact-projected -> target final error:      0.035248578750697965
Gaussian -> exact-projected final error:    0.00029022869338069174
Gaussian -> original target final error:    0.03500028070905269

target trace distance:                      0.02474893583280391
target population L2 error:                 0.02810899300694737

exact target populations:                   [0.22600611046735578, 0.7739938895326441]
Gaussian populations:                       [0.24588217003489354, 0.7541178299651065]

exact purity:                               0.6762081969769371
Gaussian purity:                            0.662382009165944
purity error:                               0.013826187810993096

exact coherence:                            [0.11414352736213967, -0.001643821331551652]
Gaussian coherence:                         [0.12888972852409641, -0.0016034206503612437]
coherence magnitude error:                  0.014744338241561974
coherence phase error / rad:                0.0019607485027196615

Bloch-vector error:                         0.049497871665607804
maximum generalized norm drift:             1.3083560634896685e-06
maximum overlap condition number:           2235.290713199147
```

## 4. Main error decomposition

For the nine-Gaussian reference,

$$
\epsilon_{\rm init}
=
0.03545458,
$$

$$
\boxed{
\epsilon_{\rm dyn}
=
0.00029022869
},
$$

and

$$
\epsilon_{\rm target}
=
0.035000281.
$$

The projected-state propagation error is over two orders of magnitude smaller than the
initial reduced-density representation error.

Therefore the compact v0.12 benchmark is currently limited predominantly by the
finite Gaussian representation of the intended coordinate-dependent initial
electronic state.

## 5. v0.11 comparison

The saved v0.11 context gives

```text
v0.11 population error:               0.012877374121210683
v0.11 center-frozen full-density error:0.15991833275047374
v0.11 trace distance:                  0.11307933752390675
v0.11 purity:                          0.643115591123172
v0.11 coherence magnitude:             0.03730993502993936
exact coherence magnitude:             0.1141553633748393
v0.11 coherence phase error / rad:     1.367544547628621
v0.11 Bloch-vector error:              0.2261586750478135
```

The v0.12 reference has a somewhat larger diagonal population error than the v0.11
optimized-spawning reference, but a much smaller full-density discrepancy and
dramatically smaller coherence phase error.

This is why the release uses the complete density matrix rather than selecting the
method with the smallest one-dimensional population error.

## 6. Acceptance

```json
{
  "checks": {
    "coherence_phase": true,
    "conditioning": true,
    "initial_density_representation": true,
    "norm": true,
    "projected_dynamics": true,
    "target_full_density": true,
    "target_population": true
  },
  "passed": true,
  "thresholds": {
    "max_coherence_phase_error": 0.01,
    "max_condition_number": 100000.0,
    "max_initial_density_error": 0.05,
    "max_norm_drift": 0.0001,
    "max_projected_dynamics_density_error": 0.001,
    "max_target_density_error": 0.05,
    "max_target_population_error": 0.05
  }
}
```

All configured v0.12 release criteria pass.

## 7. Spawning after the projected bank

An exploratory 9 -> 12 -> 15 Gaussian spawning extension was also tested during the
build.

Additional threshold-triggered spawned functions did not materially reduce the final
target-density error, while the overlap condition number increased substantially.

The release therefore does not claim that indiscriminate post-projection spawning
improves this compact benchmark.

The scientific implication is that future adaptive basis growth should be driven by a
measured representation/residual deficiency, not simply by a coupling threshold.
