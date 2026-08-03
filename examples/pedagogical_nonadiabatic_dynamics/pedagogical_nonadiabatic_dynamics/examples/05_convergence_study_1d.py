"""Verify the expected second-order global accuracy of Strang splitting."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from nonadiabatic_dynamics import (
    PeriodicGrid1D,
    build_direct_hamiltonian_1d,
    diagonalize_hamiltonian,
    diagonalize_path,
    gaussian_wavepacket_1d,
    potential_propagator_2x2,
    prepare_adiabatic_wavepacket_1d,
    propagate_from_eigendecomposition,
    smooth_single_avoided_crossing,
    split_operator_step_1d,
)
from nonadiabatic_dynamics.observables import phase_aligned_error


OUTPUT = Path(__file__).resolve().parents[1] / "outputs"


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)

    grid = PeriodicGrid1D.create(128, -14.0, 14.0)
    mass = 500.0
    final_time = 40.0
    time_steps = np.array([0.80, 0.40, 0.20, 0.10])

    potential, _ = smooth_single_avoided_crossing(grid.x)
    _, adiabatic_vectors = diagonalize_path(potential)
    scalar = gaussian_wavepacket_1d(
        grid.x,
        center=-4.0,
        width=0.75,
        momentum=4.0,
        dx=grid.dx,
    )
    psi_initial = prepare_adiabatic_wavepacket_1d(
        scalar,
        adiabatic_vectors,
        state=0,
        dx=grid.dx,
    )

    hamiltonian = build_direct_hamiltonian_1d(
        grid.k,
        mass,
        potential,
    )
    eigenvalues, eigenvectors = diagonalize_hamiltonian(hamiltonian)
    psi_reference = propagate_from_eigendecomposition(
        psi_initial,
        [final_time],
        eigenvalues,
        eigenvectors,
        spatial_shape=(grid.x.size,),
    )[0]

    fidelities = []
    errors = []

    for dt in time_steps:
        n_steps = int(round(final_time / dt))
        if not np.isclose(n_steps * dt, final_time):
            raise RuntimeError("Each time step must divide the final time.")

        potential_half = potential_propagator_2x2(
            potential,
            0.5 * dt,
        )
        kinetic = np.exp(
            -1j * grid.k**2 * dt / (2.0 * mass)
        )

        psi = psi_initial.copy()
        for _ in range(n_steps):
            psi = split_operator_step_1d(
                psi,
                potential_half,
                kinetic,
            )

        fidelity, error = phase_aligned_error(
            psi_reference,
            psi,
            grid.dx,
        )
        fidelities.append(fidelity)
        errors.append(error)

    fidelities = np.asarray(fidelities)
    errors = np.asarray(errors)

    # Fit log(error) = p log(dt) + constant.
    fitted_order = np.polyfit(
        np.log(time_steps),
        np.log(errors),
        deg=1,
    )[0]

    print("1D Strang-splitting convergence study")
    for dt, fidelity, error in zip(time_steps, fidelities, errors):
        print(
            f"dt={dt:7.4f}  fidelity={fidelity:.12f}  "
            f"phase-aligned error={error:.6e}"
        )
    print(f"Fitted global convergence order: {fitted_order:.6f}")

    figure = plt.figure(figsize=(6.5, 4.6))
    axis = figure.add_subplot(111)
    axis.loglog(time_steps, errors, marker="o", label="Measured error")

    reference_curve = errors[-1] * (time_steps / time_steps[-1]) ** 2
    axis.loglog(
        time_steps,
        reference_curve,
        linestyle="--",
        label=r"$O(\Delta t^2)$ reference",
    )
    axis.set_xlabel(r"$\Delta t$")
    axis.set_ylabel("Phase-aligned wavefunction error")
    axis.set_title("Strang split-operator time-step convergence")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "05_1d_time_step_convergence.png", dpi=180)
    plt.close(figure)

    np.savez(
        OUTPUT / "05_1d_convergence_data.npz",
        time_steps=time_steps,
        fidelities=fidelities,
        errors=errors,
        fitted_order=fitted_order,
    )


if __name__ == "__main__":
    main()
