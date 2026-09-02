from dataclasses import dataclass
import numpy as np

from .gaussian import frozen_gaussian, kinetic_on_gaussian
from .grids import inner_product
from .potentials import avoided_crossing_diabatic


@dataclass(frozen=True)
class TrajectoryBasisFunction:
    """Minimal frozen Gaussian trajectory basis function."""
    state: int
    q: float
    p: float
    alpha: float


def spawn_child(parent, target_state):
    """Create a child Gaussian on another electronic state at the same phase-space point.

    This deliberately simple choice isolates the basis-growth idea. Production
    FMS/AIMS spawning uses more sophisticated placement/optimization criteria.
    """
    if target_state == parent.state:
        raise ValueError("The child must be placed on a different electronic state.")

    return TrajectoryBasisFunction(
        state=int(target_state),
        q=float(parent.q),
        p=float(parent.p),
        alpha=float(parent.alpha),
    )


def coupling_strength_at_tbf(tbf, model=avoided_crossing_diabatic):
    H = model(float(tbf.q))
    return abs(H[0, 1])


def should_spawn(tbf, threshold=0.005, model=avoided_crossing_diabatic):
    """Simple coupling-threshold spawning criterion for the demonstration."""
    return coupling_strength_at_tbf(tbf, model=model) >= threshold


def coupled_basis_matrices(
    x,
    dx,
    basis,
    mass=1.0,
    model=avoided_crossing_diabatic,
):
    """Return S and H for a two-state diabatic Gaussian basis.

    Basis states |g_i>|I_i> are electronically orthogonal, so overlap is zero
    when I_i != I_j.

    H_ij contains kinetic energy only for equal electronic states and the
    appropriate diabatic potential matrix element V_{I_i,I_j}(x).
    """
    n = len(basis)
    S = np.zeros((n, n), dtype=complex)
    H = np.zeros((n, n), dtype=complex)

    V = model(x)

    gs = [
        frozen_gaussian(x, tbf.q, tbf.p, tbf.alpha)
        for tbf in basis
    ]

    Tgs = [
        kinetic_on_gaussian(x, tbf.q, tbf.p, tbf.alpha, mass=mass)
        for tbf in basis
    ]

    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            if bi.state == bj.state:
                S[i, j] = inner_product(gs[i], gs[j], dx)
                kinetic = inner_product(gs[i], Tgs[j], dx)
            else:
                kinetic = 0.0

            Vij = V[:, bi.state, bj.state]
            potential = inner_product(gs[i], Vij * gs[j], dx)
            H[i, j] = kinetic + potential

    return S, H


def propagate_static_basis_coefficients(C, S, H, dt):
    """One norm-stable Cayley step for a fixed nonorthogonal basis.

    i S Cdot = H C

    Crank-Nicolson:
        (S + i dt H/2) C_{n+1}
        =
        (S - i dt H/2) C_n
    """
    C = np.asarray(C, dtype=complex)

    lhs = S + 0.5j * dt * H
    rhs = (S - 0.5j * dt * H) @ C

    return np.linalg.solve(lhs, rhs)


def state_populations(C, S, basis):
    """Population by electronic state for an electronically orthogonal diabatic basis."""
    C = np.asarray(C, dtype=complex)
    populations = np.zeros(2, dtype=float)

    for state in (0, 1):
        idx = [i for i, b in enumerate(basis) if b.state == state]
        if not idx:
            continue

        block = S[np.ix_(idx, idx)]
        coeff = C[idx]
        populations[state] = np.real(np.vdot(coeff, block @ coeff))

    total = populations.sum()
    if total > 0.0:
        populations /= total

    return populations
