"""2D wavepacket dynamics near a conical intersection."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from nonadiabatic_dynamics import (
    PeriodicGrid2D,
    adiabatic_populations_2d,
    build_direct_hamiltonian_2d,
    diagonalize_hamiltonian,
    diabatic_populations,
    gaussian_wavepacket_2d,
    linear_vibronic_coupling_2d,
    lvc_analytic_eigensystem,
    norm,
    potential_propagator_2x2,
    prepare_adiabatic_wavepacket_2d,
    propagate_from_eigendecomposition,
    split_operator_step_2d,
)
from nonadiabatic_dynamics.observables import phase_aligned_error


OUTPUT = Path(__file__).resolve().parents[1] / "outputs"


def propagate_fft(
    grid: PeriodicGrid2D,
    potential: np.ndarray,
    psi_initial: np.ndarray,
    *,
    mass_x: float,
    mass_y: float,
    dt: float,
    final_time: float,
    save_every: int,
) -> tuple[np.ndarray, np.ndarray]:
    x_momentum, y_momentum = np.meshgrid(
        grid.kx,
        grid.ky,
        indexing="xy",
    )
    kinetic_energy = (
        x_momentum**2 / (2.0 * mass_x)
        + y_momentum**2 / (2.0 * mass_y)
    )
    kinetic_step = np.exp(-1j * kinetic_energy * dt)
    potential_half = potential_propagator_2x2(potential, 0.5 * dt)

    n_steps = int(round(final_time / dt))
    times = [0.0]
    wavefunctions = [psi_initial.copy()]
    psi = psi_initial.copy()

    for step in range(1, n_steps + 1):
        psi = split_operator_step_2d(
            psi,
            potential_half,
            kinetic_step,
        )
        if step % save_every == 0:
            times.append(step * dt)
            wavefunctions.append(psi.copy())

    return np.asarray(times), np.asarray(wavefunctions)


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)

    # ------------------------------------------------------------------
    # Main physically useful FFT calculation.
    # ------------------------------------------------------------------
    grid = PeriodicGrid2D.create(
        nx=64,
        ny=64,
        x_min=-8.0,
        x_max=8.0,
        y_min=-8.0,
        y_max=8.0,
    )
    x, y = np.meshgrid(grid.x, grid.y, indexing="xy")
    potential, _ = linear_vibronic_coupling_2d(x, y)
    _, adiabatic_vectors = lvc_analytic_eigensystem(x, y)

    scalar = gaussian_wavepacket_2d(
        x,
        y,
        center_x=-3.5,
        center_y=0.60,
        width_x=0.55,
        width_y=0.55,
        momentum_x=5.0,
        momentum_y=0.0,
        dx=grid.dx,
        dy=grid.dy,
    )
    psi_initial = prepare_adiabatic_wavepacket_2d(
        scalar,
        adiabatic_vectors,
        state=0,
        dx=grid.dx,
        dy=grid.dy,
    )

    times, wavefunctions = propagate_fft(
        grid,
        potential,
        psi_initial,
        mass_x=100.0,
        mass_y=100.0,
        dt=0.05,
        final_time=120.0,
        save_every=40,
    )

    norms = np.array(
        [norm(state, grid.dx * grid.dy) for state in wavefunctions]
    )
    populations_diabatic = np.array(
        [
            diabatic_populations(state, grid.dx * grid.dy)
            for state in wavefunctions
        ]
    )
    populations_adiabatic = np.array(
        [
            adiabatic_populations_2d(
                state,
                adiabatic_vectors,
                grid.dx,
                grid.dy,
            )
            for state in wavefunctions
        ]
    )

    # ------------------------------------------------------------------
    # Tiny-grid direct-diagonalization benchmark.
    #
    # This verifies the algorithm but is not a converged 2D production grid.
    # Dense diagonalization scales cubically with 2*Nx*Ny.
    # ------------------------------------------------------------------
    benchmark_grid = PeriodicGrid2D.create(
        nx=14,
        ny=14,
        x_min=-6.0,
        x_max=6.0,
        y_min=-6.0,
        y_max=6.0,
    )
    xb, yb = np.meshgrid(
        benchmark_grid.x,
        benchmark_grid.y,
        indexing="xy",
    )
    benchmark_potential, _ = linear_vibronic_coupling_2d(xb, yb)
    _, benchmark_vectors = lvc_analytic_eigensystem(xb, yb)

    benchmark_scalar = gaussian_wavepacket_2d(
        xb,
        yb,
        center_x=-2.0,
        center_y=0.70,
        width_x=0.90,
        width_y=0.90,
        momentum_x=2.0,
        momentum_y=0.0,
        dx=benchmark_grid.dx,
        dy=benchmark_grid.dy,
    )
    benchmark_initial = prepare_adiabatic_wavepacket_2d(
        benchmark_scalar,
        benchmark_vectors,
        state=0,
        dx=benchmark_grid.dx,
        dy=benchmark_grid.dy,
    )

    benchmark_times, benchmark_fft = propagate_fft(
        benchmark_grid,
        benchmark_potential,
        benchmark_initial,
        mass_x=50.0,
        mass_y=50.0,
        dt=0.02,
        final_time=8.0,
        save_every=40,
    )

    benchmark_hamiltonian = build_direct_hamiltonian_2d(
        benchmark_grid.kx,
        benchmark_grid.ky,
        50.0,
        50.0,
        benchmark_potential,
    )
    benchmark_e, benchmark_u = diagonalize_hamiltonian(
        benchmark_hamiltonian
    )
    benchmark_direct = propagate_from_eigendecomposition(
        benchmark_initial,
        benchmark_times,
        benchmark_e,
        benchmark_u,
        spatial_shape=(benchmark_grid.y.size, benchmark_grid.x.size),
    )

    fidelity, aligned_error = phase_aligned_error(
        benchmark_direct[-1],
        benchmark_fft[-1],
        benchmark_grid.dx * benchmark_grid.dy,
    )

    print("2D conical-intersection wavepacket dynamics")
    print(f"Main FFT norm error: {np.max(np.abs(norms - 1.0)):.6e}")
    print("Final FFT diabatic populations:", populations_diabatic[-1])
    print("Final FFT adiabatic populations:", populations_adiabatic[-1])
    print(f"Tiny-grid direct benchmark fidelity: {fidelity:.12f}")
    print(f"Tiny-grid phase-aligned L2 error: {aligned_error:.6e}")

    density_initial = np.sum(np.abs(wavefunctions[0]) ** 2, axis=0)
    density_middle = np.sum(
        np.abs(wavefunctions[len(wavefunctions) // 2]) ** 2,
        axis=0,
    )
    density_final = np.sum(np.abs(wavefunctions[-1]) ** 2, axis=0)

    for label, density in [
        ("initial", density_initial),
        ("middle", density_middle),
        ("final", density_final),
    ]:
        figure = plt.figure(figsize=(6.2, 5.2))
        axis = figure.add_subplot(111)
        contour = axis.contourf(x, y, density, levels=24)
        figure.colorbar(contour, ax=axis, label=r"$|\Psi(x,y,t)|^2$")
        axis.set_xlabel(r"$x$")
        axis.set_ylabel(r"$y$")
        axis.set_aspect("equal")
        axis.set_title(f"2D nuclear density: {label}")
        figure.tight_layout()
        figure.savefig(OUTPUT / f"04_2d_density_{label}.png", dpi=180)
        plt.close(figure)

    figure = plt.figure(figsize=(7.4, 4.6))
    axis = figure.add_subplot(111)
    axis.plot(times, populations_adiabatic[:, 0], label="Lower adiabatic")
    axis.plot(times, populations_adiabatic[:, 1], label="Upper adiabatic")
    axis.set_xlabel("Time")
    axis.set_ylabel("Population")
    axis.set_title("2D nonadiabatic population transfer")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "04_2d_adiabatic_populations.png", dpi=180)
    plt.close(figure)

    np.savez(
        OUTPUT / "04_2d_dynamics_data.npz",
        x=grid.x,
        y=grid.y,
        times=times,
        wavefunctions=wavefunctions,
        norms=norms,
        diabatic_populations=populations_diabatic,
        adiabatic_populations=populations_adiabatic,
        benchmark_times=benchmark_times,
        benchmark_fft=benchmark_fft,
        benchmark_direct=benchmark_direct,
        benchmark_fidelity=fidelity,
        benchmark_phase_aligned_error=aligned_error,
    )


if __name__ == "__main__":
    main()
