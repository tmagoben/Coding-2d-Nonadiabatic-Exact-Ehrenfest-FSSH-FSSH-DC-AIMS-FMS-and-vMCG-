"""Regression tests for the pedagogical algorithms."""

import numpy as np

from nonadiabatic_dynamics import (
    PeriodicGrid1D,
    PeriodicGrid2D,
    build_direct_hamiltonian_1d,
    build_direct_hamiltonian_2d,
    diagonalize_hamiltonian,
    diagonalize_path,
    gaussian_wavepacket_1d,
    gaussian_wavepacket_2d,
    linear_vibronic_coupling_2d,
    lvc_analytic_eigensystem,
    norm,
    potential_propagator_2x2,
    prepare_adiabatic_wavepacket_1d,
    prepare_adiabatic_wavepacket_2d,
    propagate_from_eigendecomposition,
    split_operator_step_1d,
    split_operator_step_2d,
    smooth_single_avoided_crossing,
)
from nonadiabatic_dynamics.observables import phase_aligned_error


def test_local_potential_propagator_is_unitary() -> None:
    x = np.linspace(-2.0, 2.0, 31)
    potential, _ = smooth_single_avoided_crossing(x)
    propagator = potential_propagator_2x2(potential, 0.37)
    identity = np.einsum(
        "xba,xbc->xac",
        propagator.conj(),
        propagator,
    )
    assert np.allclose(identity, np.eye(2)[None], atol=1.0e-13)


def test_1d_split_operator_matches_direct_diagonalization() -> None:
    grid = PeriodicGrid1D.create(64, -12.0, 12.0)
    mass = 300.0
    dt = 0.02
    final_time = 1.0

    potential, _ = smooth_single_avoided_crossing(grid.x)
    _, vectors = diagonalize_path(potential)
    scalar = gaussian_wavepacket_1d(
        grid.x,
        center=-3.0,
        width=0.70,
        momentum=3.0,
        dx=grid.dx,
    )
    psi0 = prepare_adiabatic_wavepacket_1d(
        scalar,
        vectors,
        state=0,
        dx=grid.dx,
    )

    potential_half = potential_propagator_2x2(potential, 0.5 * dt)
    kinetic = np.exp(-1j * grid.k**2 * dt / (2.0 * mass))
    psi_fft = psi0.copy()
    for _ in range(int(round(final_time / dt))):
        psi_fft = split_operator_step_1d(
            psi_fft,
            potential_half,
            kinetic,
        )

    hamiltonian = build_direct_hamiltonian_1d(grid.k, mass, potential)
    energies, eigenvectors = diagonalize_hamiltonian(hamiltonian)
    psi_direct = propagate_from_eigendecomposition(
        psi0,
        [final_time],
        energies,
        eigenvectors,
        spatial_shape=(grid.x.size,),
    )[0]

    fidelity, error = phase_aligned_error(
        psi_direct,
        psi_fft,
        grid.dx,
    )
    assert fidelity > 1.0 - 1.0e-11
    assert error < 5.0e-6
    assert abs(norm(psi_fft, grid.dx) - 1.0) < 1.0e-12


def test_2d_split_operator_matches_tiny_direct_diagonalization() -> None:
    grid = PeriodicGrid2D.create(8, 8, -4.0, 4.0, -4.0, 4.0)
    x, y = np.meshgrid(grid.x, grid.y, indexing="xy")
    potential, _ = linear_vibronic_coupling_2d(x, y)
    _, vectors = lvc_analytic_eigensystem(x, y)

    scalar = gaussian_wavepacket_2d(
        x,
        y,
        center_x=-1.0,
        center_y=0.5,
        width_x=0.8,
        width_y=0.8,
        momentum_x=1.0,
        momentum_y=0.0,
        dx=grid.dx,
        dy=grid.dy,
    )
    psi0 = prepare_adiabatic_wavepacket_2d(
        scalar,
        vectors,
        state=0,
        dx=grid.dx,
        dy=grid.dy,
    )

    mass = 40.0
    dt = 0.01
    final_time = 0.20
    potential_half = potential_propagator_2x2(potential, 0.5 * dt)
    kx, ky = np.meshgrid(grid.kx, grid.ky, indexing="xy")
    kinetic = np.exp(
        -1j * (kx**2 + ky**2) * dt / (2.0 * mass)
    )

    psi_fft = psi0.copy()
    for _ in range(int(round(final_time / dt))):
        psi_fft = split_operator_step_2d(
            psi_fft,
            potential_half,
            kinetic,
        )

    hamiltonian = build_direct_hamiltonian_2d(
        grid.kx,
        grid.ky,
        mass,
        mass,
        potential,
    )
    energies, eigenvectors = diagonalize_hamiltonian(hamiltonian)
    psi_direct = propagate_from_eigendecomposition(
        psi0,
        [final_time],
        energies,
        eigenvectors,
        spatial_shape=(grid.y.size, grid.x.size),
    )[0]

    fidelity, error = phase_aligned_error(
        psi_direct,
        psi_fft,
        grid.dx * grid.dy,
    )
    assert fidelity > 1.0 - 1.0e-11
    assert error < 5.0e-6
