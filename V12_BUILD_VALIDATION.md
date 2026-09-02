# v0.12 Build Validation Report

Validated on 2026-08-12.

## Source validation

All 182 Python files in the repository were parsed successfully with Python's AST
parser before packaging.

## Automated regression suite

```text
149 passed in 6.73 s
```

The cumulative suite retains all v0.1-v0.11 tests and adds v0.12 checks for:

- exact analytic unequal-width LVC Gaussian potential and total Hamiltonian matrix
  elements against direct 2D quadrature;
- Hermiticity of exact LVC Gaussian S/H matrices;
- center adiabatic spinor derivatives against finite differences;
- the continuous moving-basis identity `Sdot = T + T^dagger`;
- generalized midpoint Cayley reduction to fixed-basis Cayley propagation;
- coherence magnitude, wrapped phase, trace-distance, and Bloch-vector diagnostics;
- explicit local-diabatic spinor transport;
- pair-preserving pruning of a spinor-complete nuclear Gaussian basis;
- spinor-complete exact-LVC S/H blocks and short-time norm conservation;
- coordinate-dependent Born-Huang TBF reconstruction of the exact initial packet;
- Born-Huang projected S/H Hermiticity;
- Born-Huang moving-basis finite-difference consistency;
- short Born-Huang projected propagation;
- exact-grid initial-state projection into a spinor-complete Gaussian bank;
- independent v0.12 acceptance checks for representation and dynamics errors.

## Representative executable examples

### Exact LVC Gaussian matrices

```text
S =
[[1.        +0.j         0.50728213-0.07057896j]
 [0.50728213+0.07057896j 1.        +0.j        ]]

H =
[[0.11207225+0.j        0.02370716-0.0015355j]
 [0.02370716+0.0015355j 0.0696658 +0.j       ]]

||S-S^dag|| = 0.0
||H-H^dag|| = 4.9065389333867974e-18
cond(S) = 3.099776095083131
```

### Initial representation audit

```text
exact coordinate-dependent initial populations:
[0.2258137522757219, 0.774186247724278]

exact coordinate-dependent initial purity:
0.6764597760317345

center-frozen initial populations:
[0.03846153846153845, 0.9615384615384617]

center-frozen initial purity:
0.9999999999999998

center-frozen reduced-density error:
0.28703562527170995
```

This confirms that the old center-frozen and exact coordinate-dependent initial
ansätze are already different at t=0.

### Nine-Gaussian representation-consistent benchmark

```text
initial projection fidelity:           0.832276023595292
initial reduced-density error:         0.03545457994295867
projected-state dynamics error:        0.00029022869338069174
original-target final density error:   0.03500028070905269
population L2 error:                   0.02810899300694737
trace distance:                        0.02474893583280391
purity error:                          0.013826187810993096
coherence phase error / rad:           0.0019607485027196615
maximum generalized norm drift:        1.3083560634896685e-06
maximum overlap condition number:      2235.290713199147
```

### Born-Huang initial-state reconstruction

```text
wavefunction L2 difference: 0.0
```

One coordinate-dependent `g(R) Phi_a(R)` TBF therefore reproduces the intended
localized adiabatic initial packet exactly on the same discrete benchmark grid.

## v0.11 -> v0.12 density/coherence comparison

```text
v0.11 full-density error:          0.15991833275047374
v0.12 full-density error:          0.03500028070905269

v0.11 coherence phase error / rad: 1.367544547628621
v0.12 coherence phase error / rad: 0.0019607485027196615

v0.11 population error:            0.012877374121210683
v0.12 population error:            0.02810899300694737
```

v0.12 is not selected merely because one population number is smaller.  Its release
criterion is the representation-consistent complete reduced electronic density.

## Acceptance result

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

All configured release criteria pass.

## Scientific interpretation

The key v0.12 number is

```text
projected-state dynamics density error =
0.00029022869338069174
```

because the Gaussian and exact calculations compared there begin from the identical
projected initial wavefunction.

The original-target error

```text
0.03500028070905269
```

is much larger and closely follows the finite initial-state representation error

```text
0.03545457994295867.
```

Therefore the compact release benchmark identifies initial representation of the
coordinate-dependent electronic state as the leading remaining approximation.

## PySCF status

The inherited PySCF SA-CASSCF backend, many-electron overlap tracking, and gauge-graph
layers remain included.

The v0.12 exact-LVC and Born-Huang release benchmarks are analytic-model
implementations.  No claim is made that an exact molecular/global-diabatic PySCF
Hamiltonian has been implemented.

See `V12_PYSCF_BRIDGE.md`.
