# Validation report

Validated in the build environment on 2026-08-12.

## Automated tests

```text
13 passed
```

The tests cover:

- Gaussian normalization and moments;
- analytic versus numerical overlap;
- split-operator norm conservation;
- harmonic Heller center dynamics;
- Heller TGA versus exact harmonic propagation;
- moving-basis Hermiticity;
- the identity dS/dt = tau^dagger + tau;
- moving-basis norm conservation;
- McLachlan TDVP residual reduction;
- spawning basis construction;
- fixed-basis Cayley norm conservation.

## Representative numerical checks

Harmonic Heller example:

```text
exact norm      = 1.000000000000005
TGA norm        = 1.0
fidelity        = 0.9999999999994398
```

Anharmonic quartic example:

```text
exact vs TGA fidelity = 0.9863215243196044
```

This lower fidelity is expected because the exact anharmonic wavepacket develops
non-Gaussian structure.

Moving two-Gaussian basis:

```text
initial norm = 1.0000000000000002
final norm   = 1.0000000000000064
max |N-1|    = 1.94e-14
```

Variational two-Gaussian tangent projection:

```text
projected TDSE residual = 0.0171030
zero-velocity residual  = 1.6738276
```

The purpose of these values is regression/consistency checking; they are not presented
as molecular benchmark data.
