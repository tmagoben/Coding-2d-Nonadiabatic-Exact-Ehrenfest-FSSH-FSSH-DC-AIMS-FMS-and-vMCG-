# How this repository connects to the companion nonadiabatic-dynamics projects

The intended research/software progression is:

```text
Exact grid dynamics
        |
        +--> Ehrenfest
        |
        +--> FSSH
        |
        +--> FSSH + decoherence
        |
        v
Gaussian nuclear wavefunctions
        |
        +--> Heller single Gaussian
        |
        +--> moving Gaussian basis
        |
        +--> variational multi-Gaussian TDVP
        |
        +--> FMS/AIMS spawning foundation
```

The key conceptual transition is from

```text
classical nuclear trajectories
```

to

```text
time-dependent nuclear basis functions carrying complex amplitudes.
```

A Gaussian center may still follow a trajectory, but the physical object is the
coherent wavefunction

$$
\Psi=\sum_j C_j g_j,
$$

not the collection of centers by itself.

## Recommended future extension order

1. two-state exact split-operator benchmark;
2. adiabatic multi-state frozen Gaussian basis;
3. explicit derivative-coupling matrix elements;
4. gauge/phase-consistent overlap transport;
5. dynamic FMS spawning with parent/child basis functions;
6. AIMS electronic-structure interface;
7. 2D conical-intersection model;
8. connection to GeneralDIA/PySCF reference data.

This order keeps each new scientific complication isolated and testable.
