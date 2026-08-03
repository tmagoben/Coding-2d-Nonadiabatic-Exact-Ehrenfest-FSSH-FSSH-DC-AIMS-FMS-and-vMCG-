# Pedagogical Nonadiabatic Quantum Dynamics

This directory develops a two-state nonadiabatic problem sequentially:

1. construct diabatic model Hamiltonians
2. diagonalize them into an adiabatic representation
3. track electronic states and phases
4. calculate derivative couplings
5. integrate a pathwise adiabatic-to-diabatic transformation
6. prepare nuclear Gaussian wavepackets
7. propagate coupled dynamics with an FFT split operator
8. construct the complete finite-grid Hamiltonian
9. propagate the same problem by direct diagonalization
10. compare norms, wavefunctions, and electronic populations
11. repeat the construction in two nuclear dimensions around a conical
    intersection.

Read [`DYNAMICS_TUTORIAL.md`](DYNAMICS_TUTORIAL.md) for the derivation, then [`IMPLEMENTATION_WALKTHROUGH.md`](IMPLEMENTATION_WALKTHROUGH.md) for the line-by-line mapping into Python.

## Install

```bash
python -m pip install -r requirements.txt
python -m pip install -e . --no-build-isolation
```

## Run sequentially

```bash
python examples/01_static_1d_derivative_coupling_and_adt.py
python examples/02_dynamics_1d_fft_vs_diagonalization.py
python examples/03_static_2d_conical_intersection.py
python examples/04_dynamics_2d_fft_and_direct_benchmark.py
python examples/05_convergence_study_1d.py
pytest
```

Generated figures and arrays are written to `outputs/`.

## Why use both propagators?

The FFT split operator scales efficiently to large grids but introduces
time-step splitting error. Full diagonalization is exact for the selected
finite grid but scales cubically in the total basis dimension(AKA the "curse of dimensionality). It is therefore
used as a small-system reference implementation. 
