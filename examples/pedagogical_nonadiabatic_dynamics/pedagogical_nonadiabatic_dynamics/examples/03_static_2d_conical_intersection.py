"""Static 2D example: conical-intersection NAC field and Berry phase."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from nonadiabatic_dynamics import (
    hellmann_feynman_derivative_couplings,
    integrate_adt_path,
    linear_vibronic_coupling_2d,
    lvc_analytic_derivative_coupling,
    lvc_analytic_eigensystem,
)


OUTPUT = Path(__file__).resolve().parents[1] / "outputs"


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)

    coordinate = np.linspace(-3.0, 3.0, 121)
    x, y = np.meshgrid(coordinate, coordinate, indexing="xy")
    energies, vectors = lvc_analytic_eigensystem(x, y)
    _, derivatives = linear_vibronic_coupling_2d(x, y)

    tau = hellmann_feynman_derivative_couplings(
        energies.reshape(-1, 2),
        vectors.reshape(-1, 2, 2),
        derivatives.reshape(-1, 2, 2, 2),
        gap_tolerance=1.0e-10,
    ).reshape(x.shape + (2, 2, 2))

    exact_x, exact_y, valid = lvc_analytic_derivative_coupling(
        x,
        y,
        singular_radius=1.0e-8,
    )
    error = max(
        np.nanmax(np.abs(tau[..., 0, 0, 1][valid] - exact_x[valid])),
        np.nanmax(np.abs(tau[..., 1, 0, 1][valid] - exact_y[valid])),
    )

    phi = np.linspace(0.0, 2.0 * np.pi, 2401)
    radius = 1.5
    x_loop = radius * np.cos(phi)
    y_loop = radius * np.sin(phi)
    e_loop, u_loop = lvc_analytic_eigensystem(x_loop, y_loop)
    _, dh_loop = linear_vibronic_coupling_2d(x_loop, y_loop)

    dx_dphi = -radius * np.sin(phi)
    dy_dphi = radius * np.cos(phi)
    dh_dphi = (
        dh_loop[:, 0] * dx_dphi[:, None, None]
        + dh_loop[:, 1] * dy_dphi[:, None, None]
    )
    tau_phi = hellmann_feynman_derivative_couplings(
        e_loop,
        u_loop,
        dh_dphi[:, None],
    )[:, 0]

    adt = integrate_adt_path(phi, tau_phi, u_loop[0].T)
    holonomy = adt[-1] @ adt[0].T
    berry_integral = trapezoid(tau_phi[:, 0, 1], phi)

    print("Static 2D conical intersection")
    print(f"Analytic-versus-HF NAC error: {error:.6e}")
    print(f"Closed-loop NAC integral: {berry_integral:.12f}")
    print("Holonomy:")
    print(holonomy)

    gap = energies[..., 1] - energies[..., 0]
    stride = 8

    figure = plt.figure(figsize=(7.0, 5.6))
    axis = figure.add_subplot(111)
    contour = axis.contourf(x, y, gap, levels=24)
    figure.colorbar(contour, ax=axis, label=r"$E_1-E_0$")
    axis.quiver(
        x[::stride, ::stride],
        y[::stride, ::stride],
        tau[::stride, ::stride, 0, 0, 1],
        tau[::stride, ::stride, 1, 0, 1],
    )
    axis.set_xlabel(r"$x$")
    axis.set_ylabel(r"$y$")
    axis.set_aspect("equal")
    axis.set_title("Conical-intersection gap and NAC vector field")
    figure.tight_layout()
    figure.savefig(OUTPUT / "03_2d_gap_and_nac.png", dpi=180)
    plt.close(figure)

    np.savez(
        OUTPUT / "03_static_2d_data.npz",
        x=x,
        y=y,
        energies=energies,
        derivative_couplings=tau,
        phi=phi,
        tau_phi=tau_phi,
        adt=adt,
        holonomy=holonomy,
        berry_integral=berry_integral,
    )


if __name__ == "__main__":
    main()
