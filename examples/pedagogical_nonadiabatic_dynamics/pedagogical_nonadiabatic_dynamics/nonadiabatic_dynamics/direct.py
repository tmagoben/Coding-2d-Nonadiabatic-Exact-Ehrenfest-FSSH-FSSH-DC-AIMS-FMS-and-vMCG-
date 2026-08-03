"""Direct finite-grid Hamiltonian construction and exact propagation."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def _spectral_kinetic_matrix_1d(
    momenta: FloatArray,
    mass: float,
) -> ComplexArray:
    """Return the periodic Fourier-spectral kinetic matrix."""
    if mass <= 0.0:
        raise ValueError("mass must be positive.")
    kinetic_symbol = momenta**2 / (2.0 * mass)
    identity = np.eye(momenta.size, dtype=complex)
    matrix = np.fft.ifft(
        kinetic_symbol[:, None] * np.fft.fft(identity, axis=0),
        axis=0,
    )
    return 0.5 * (matrix + matrix.conj().T)


def build_direct_hamiltonian_1d(
    momenta: ArrayLike,
    mass: float,
    potential: ArrayLike,
) -> ComplexArray:
    r"""
    Build the full two-state periodic-grid Hamiltonian.

    The flattened ordering is state-major:

        [Psi_0(x_0...x_N-1), Psi_1(x_0...x_N-1)].
    """
    k = np.asarray(momenta, dtype=float)
    v = np.asarray(potential, dtype=float)
    if v.shape != (k.size, 2, 2):
        raise ValueError("Potential must have shape (n_grid,2,2).")

    n = k.size
    kinetic = _spectral_kinetic_matrix_1d(k, mass)
    hamiltonian = np.kron(np.eye(2), kinetic)

    indices = np.arange(n)
    for a in range(2):
        for b in range(2):
            hamiltonian[a * n + indices, b * n + indices] += v[:, a, b]

    return 0.5 * (hamiltonian + hamiltonian.conj().T)


def build_direct_hamiltonian_2d(
    kx: ArrayLike,
    ky: ArrayLike,
    mass_x: float,
    mass_y: float,
    potential: ArrayLike,
) -> ComplexArray:
    r"""
    Build the full two-state 2D periodic-grid Hamiltonian.

    Spatial flattening uses C order, so x is the fastest index.
    """
    kx = np.asarray(kx, dtype=float)
    ky = np.asarray(ky, dtype=float)
    v = np.asarray(potential, dtype=float)

    ny, nx = v.shape[:2]
    if v.shape != (ky.size, kx.size, 2, 2):
        raise ValueError("Potential must have shape (ny,nx,2,2).")
    if mass_x <= 0.0 or mass_y <= 0.0:
        raise ValueError("Masses must be positive.")

    tx = _spectral_kinetic_matrix_1d(kx, mass_x)
    ty = _spectral_kinetic_matrix_1d(ky, mass_y)
    spatial_kinetic = (
        np.kron(np.eye(ny), tx)
        + np.kron(ty, np.eye(nx))
    )

    n_space = nx * ny
    hamiltonian = np.kron(np.eye(2), spatial_kinetic)
    flat_v = v.reshape(n_space, 2, 2)
    indices = np.arange(n_space)

    for a in range(2):
        for b in range(2):
            hamiltonian[
                a * n_space + indices,
                b * n_space + indices,
            ] += flat_v[:, a, b]

    return 0.5 * (hamiltonian + hamiltonian.conj().T)


def diagonalize_hamiltonian(
    hamiltonian: ArrayLike,
) -> tuple[FloatArray, ComplexArray]:
    """Return all eigenvalues and eigenvectors of the finite Hamiltonian."""
    h = np.asarray(hamiltonian, dtype=complex)
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise ValueError("Hamiltonian must be square.")
    energies, vectors = np.linalg.eigh(h)
    return energies.real, vectors


def propagate_from_eigendecomposition(
    psi_initial: ArrayLike,
    times: ArrayLike,
    eigenvalues: ArrayLike,
    eigenvectors: ArrayLike,
    *,
    spatial_shape: tuple[int, ...],
) -> ComplexArray:
    r"""
    Apply the exact finite-basis propagator

        Psi(t) = W exp(-i epsilon t) W^\dagger Psi(0).

    Returns shape ``(n_times,2)+spatial_shape``.
    """
    psi0 = np.asarray(psi_initial, dtype=complex)
    time_array = np.asarray(times, dtype=float)
    energies = np.asarray(eigenvalues, dtype=float)
    vectors = np.asarray(eigenvectors, dtype=complex)

    expected_shape = (2,) + tuple(spatial_shape)
    if psi0.shape != expected_shape:
        raise ValueError(
            f"Expected initial wavefunction shape {expected_shape}, "
            f"received {psi0.shape}."
        )

    flattened = psi0.reshape(-1)
    coefficients = vectors.conj().T @ flattened
    result = np.empty((time_array.size,) + expected_shape, dtype=complex)

    for index, time in enumerate(time_array):
        propagated = vectors @ (
            np.exp(-1j * energies * time) * coefficients
        )
        result[index] = propagated.reshape(expected_shape)
    return result
