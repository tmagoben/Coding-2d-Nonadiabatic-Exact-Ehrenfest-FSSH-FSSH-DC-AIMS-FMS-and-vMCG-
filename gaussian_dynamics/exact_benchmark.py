import numpy as np

from .ci2d import diabatic_potential_2d, analytic_adiabatic_vectors
from .exact2d import run_exact_2d
from .gaussian_nd import gaussian_nd


def localized_adiabatic_packet_2d(points, q0, p0, A, state=1):
    """Construct g(R) phi_state(R) in the diabatic representation.

    The analytic real adiabatic frame carries the expected branch/gauge structure.
    For a packet localized away from the exact CI this is a transparent benchmark
    initialization.  The exact CI point must be excluded from the grid.
    """
    points = np.asarray(points, dtype=float)
    g = gaussian_nd(points, q0, p0, A)
    psi = np.zeros(points.shape[:-1] + (2,), dtype=complex)

    flat = points.reshape(-1, 2)
    out = psi.reshape(-1, 2)
    gf = g.reshape(-1)

    for i, R in enumerate(flat):
        U = analytic_adiabatic_vectors(R)
        out[i] = gf[i] * U[:, state]

    return psi


def adiabatic_populations_from_diabatic(psi, points, dx, dy):
    psi = np.asarray(psi, dtype=complex)
    points = np.asarray(points, dtype=float)
    if psi.shape[:-1] != points.shape[:-1] or psi.shape[-1] != 2:
        raise ValueError("psi and points have incompatible shapes")

    pops = np.zeros(2, dtype=float)
    flat_R = points.reshape(-1, 2)
    flat_psi = psi.reshape(-1, 2)

    for R, vector in zip(flat_R, flat_psi):
        U = analytic_adiabatic_vectors(R)
        c = U.conj().T @ vector
        pops += np.abs(c) ** 2 * dx * dy

    return pops


def run_exact_ci_reference(
    q0=np.array([0.55, 0.45]),
    p0=np.array([0.6, 0.8]),
    A=1.2*np.eye(2),
    state=1,
    mass=20.0,
    grid_n=48,
    half_width=4.0,
    dt=0.001,
    final_time=0.01,
):
    n = int(grid_n)
    dx = 2.0 * half_width / n
    x = -half_width + (np.arange(n) + 0.5) * dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    points = np.stack([X, Y], axis=-1)

    psi0 = localized_adiabatic_packet_2d(points, q0, p0, A, state=state)
    V = diabatic_potential_2d(X, Y)
    steps = int(round(final_time / dt))

    out = run_exact_2d(
        psi0,
        dx,
        dx,
        V,
        mass=mass,
        dt=dt,
        steps=steps,
        store_every=steps,
    )
    final_pops = adiabatic_populations_from_diabatic(
        out["psi"][-1], points, dx, dx
    )

    return {
        "final_populations_adiabatic": final_pops,
        "norm": float(out["norm"][-1]),
        "psi_final": out["psi"][-1],
        "points": points,
        "dx": dx,
        "dt": dt,
        "steps": steps,
    }
