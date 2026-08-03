"""1D wavepacket dynamics: FFT split operator versus direct diagonalization."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from nonadiabatic_dynamics import (
    PeriodicGrid1D,
    adiabatic_populations_1d,
    build_direct_hamiltonian_1d,
    diagonalize_hamiltonian,
    diagonalize_path,
    diabatic_populations,
    gaussian_wavepacket_1d,
    norm,
    potential_propagator_2x2,
    prepare_adiabatic_wavepacket_1d,
    propagate_from_eigendecomposition,
    split_operator_step_1d,
    smooth_single_avoided_crossing,
)
from nonadiabatic_dynamics.observables import phase_aligned_error


OUTPUT = Path(__file__).resolve().parents[1] / "outputs"


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)

    # A moderate grid is deliberately used so complete diagonalization remains
    # feasible while the FFT method is already representative.
    grid = PeriodicGrid1D.create(256, -20.0, 20.0)
    mass = 1000.0
    dt = 0.20
    final_time = 1000.0
    save_every = 50
    n_steps = int(round(final_time / dt))

    potential, _ = smooth_single_avoided_crossing(grid.x)
    energies_ad, vectors_ad = diagonalize_path(potential)

    scalar = gaussian_wavepacket_1d(
        grid.x,
        center=-8.0,
        width=0.70,
        momentum=12.0,
        dx=grid.dx,
    )
    psi_initial = prepare_adiabatic_wavepacket_1d(
        scalar,
        vectors_ad,
        state=0,
        dx=grid.dx,
    )

    potential_half = potential_propagator_2x2(potential, 0.5 * dt)
    kinetic_full = np.exp(
        -1j * grid.k**2 * dt / (2.0 * mass)
    )

    saved_times = [0.0]
    saved_fft = [psi_initial.copy()]
    psi = psi_initial.copy()

    for step in range(1, n_steps + 1):
        psi = split_operator_step_1d(
            psi,
            potential_half,
            kinetic_full,
        )
        if step % save_every == 0:
            saved_times.append(step * dt)
            saved_fft.append(psi.copy())

    saved_times = np.asarray(saved_times)
    saved_fft = np.asarray(saved_fft)

    # Direct finite-grid solution: exp(-iHt) constructed from all eigenpairs.
    hamiltonian = build_direct_hamiltonian_1d(
        grid.k,
        mass,
        potential,
    )
    eigenvalues, eigenvectors = diagonalize_hamiltonian(hamiltonian)
    saved_direct = propagate_from_eigendecomposition(
        psi_initial,
        saved_times,
        eigenvalues,
        eigenvectors,
        spatial_shape=(grid.x.size,),
    )

    fft_norms = np.array([norm(state, grid.dx) for state in saved_fft])
    direct_norms = np.array([norm(state, grid.dx) for state in saved_direct])

    fft_diabatic = np.array(
        [diabatic_populations(state, grid.dx) for state in saved_fft]
    )
    fft_adiabatic = np.array(
        [
            adiabatic_populations_1d(state, vectors_ad, grid.dx)
            for state in saved_fft
        ]
    )
    direct_adiabatic = np.array(
        [
            adiabatic_populations_1d(state, vectors_ad, grid.dx)
            for state in saved_direct
        ]
    )

    fidelity, aligned_error = phase_aligned_error(
        saved_direct[-1],
        saved_fft[-1],
        grid.dx,
    )

    print("1D FFT versus direct diagonalization")
    print(f"Grid points: {grid.x.size}")
    print(f"Final time: {final_time:.3f}")
    print(f"FFT norm error: {np.max(np.abs(fft_norms - 1.0)):.6e}")
    print(f"Direct norm error: {np.max(np.abs(direct_norms - 1.0)):.6e}")
    print(f"Final-state fidelity: {fidelity:.12f}")
    print(f"Final phase-aligned L2 error: {aligned_error:.6e}")
    print("Final FFT adiabatic populations:", fft_adiabatic[-1])
    print("Final direct adiabatic populations:", direct_adiabatic[-1])

    density_initial = np.sum(np.abs(saved_fft[0]) ** 2, axis=0)
    density_middle = np.sum(
        np.abs(saved_fft[len(saved_fft) // 2]) ** 2,
        axis=0,
    )
    density_final = np.sum(np.abs(saved_fft[-1]) ** 2, axis=0)

    figure = plt.figure(figsize=(7.4, 4.6))
    axis = figure.add_subplot(111)
    axis.plot(grid.x, density_initial, label=r"$t=0$")
    axis.plot(
        grid.x,
        density_middle,
        label=rf"$t={saved_times[len(saved_times)//2]:.0f}$",
    )
    axis.plot(grid.x, density_final, label=rf"$t={final_time:.0f}$")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$|\Psi(x,t)|^2$")
    axis.set_title("1D nuclear probability density")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "02_1d_density_snapshots.png", dpi=180)
    plt.close(figure)

    figure = plt.figure(figsize=(7.4, 4.6))
    axis = figure.add_subplot(111)
    axis.plot(saved_times, fft_adiabatic[:, 0], label="FFT lower")
    axis.plot(saved_times, fft_adiabatic[:, 1], label="FFT upper")
    axis.plot(
        saved_times,
        direct_adiabatic[:, 0],
        linestyle="--",
        label="Direct lower",
    )
    axis.plot(
        saved_times,
        direct_adiabatic[:, 1],
        linestyle="--",
        label="Direct upper",
    )
    axis.set_xlabel("Time")
    axis.set_ylabel("Adiabatic population")
    axis.set_title("Nonadiabatic population transfer")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "02_1d_populations.png", dpi=180)
    plt.close(figure)

    np.savez(
        OUTPUT / "02_1d_dynamics_data.npz",
        x=grid.x,
        times=saved_times,
        fft_wavefunctions=saved_fft,
        direct_wavefunctions=saved_direct,
        fft_diabatic_populations=fft_diabatic,
        fft_adiabatic_populations=fft_adiabatic,
        direct_adiabatic_populations=direct_adiabatic,
        fft_norms=fft_norms,
        direct_norms=direct_norms,
        fidelity=fidelity,
        phase_aligned_error=aligned_error,
    )


if __name__ == "__main__":
    main()
