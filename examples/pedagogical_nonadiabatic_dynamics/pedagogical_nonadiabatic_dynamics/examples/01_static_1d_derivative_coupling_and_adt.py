"""Static 1D example: states, derivative couplings, and diabatization."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from nonadiabatic_dynamics import (
    diagonalize_path,
    finite_difference_derivative_couplings_1d,
    hellmann_feynman_derivative_couplings,
    integrate_adt_path,
    transform_adiabatic_to_diabatic,
    smooth_single_avoided_crossing,
)


OUTPUT = Path(__file__).resolve().parents[1] / "outputs"


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)

    x = np.linspace(-8.0, 8.0, 6401)
    original_diabatic, derivative = smooth_single_avoided_crossing(x)
    energies, vectors = diagonalize_path(original_diabatic)

    tau_hf = hellmann_feynman_derivative_couplings(
        energies,
        vectors,
        derivative,
    )[:, 0]
    tau_fd = finite_difference_derivative_couplings_1d(x, vectors)

    adt = integrate_adt_path(x, tau_hf, vectors[0].T)
    recovered_diabatic = transform_adiabatic_to_diabatic(energies, adt)

    fd_error = np.max(
        np.abs(tau_hf[2:-2, 0, 1] - tau_fd[2:-2, 0, 1])
    )
    reconstruction_error = np.max(
        np.abs(recovered_diabatic - original_diabatic)
    )

    print("Static 1D avoided crossing")
    print(f"HF-versus-finite-difference NAC error: {fd_error:.6e}")
    print(f"Diabatic reconstruction error: {reconstruction_error:.6e}")

    figure = plt.figure(figsize=(7.2, 4.5))
    axis = figure.add_subplot(111)
    axis.plot(x, energies[:, 0], label=r"$E_0$")
    axis.plot(x, energies[:, 1], label=r"$E_1$")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel("Energy")
    axis.set_title("Adiabatic energies")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "01_static_1d_energies.png", dpi=180)
    plt.close(figure)

    figure = plt.figure(figsize=(7.2, 4.5))
    axis = figure.add_subplot(111)
    axis.plot(x, tau_hf[:, 0, 1], label="Hellmann-Feynman")
    axis.plot(x, tau_fd[:, 0, 1], linestyle="--", label="Finite difference")
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$\tau_{01}(x)$")
    axis.set_title("Derivative coupling")
    axis.legend()
    figure.tight_layout()
    figure.savefig(OUTPUT / "01_static_1d_nac.png", dpi=180)
    plt.close(figure)

    np.savez(
        OUTPUT / "01_static_1d_data.npz",
        x=x,
        original_diabatic=original_diabatic,
        adiabatic_energies=energies,
        adiabatic_vectors=vectors,
        tau_hf=tau_hf,
        tau_fd=tau_fd,
        adt=adt,
        recovered_diabatic=recovered_diabatic,
    )


if __name__ == "__main__":
    main()
