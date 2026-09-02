from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class LVC2DParameters:
    """Two-state linear-vibronic-coupling conical-intersection model."""
    kappa: float = 0.02
    lam: float = 0.02
    omega: float = 0.02


def diabatic_potential_2d(x, y, p=LVC2DParameters()):
    """Return V_d(x,y) with final matrix axes (...,2,2)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x, y = np.broadcast_arrays(x, y)

    common = 0.5 * p.omega**2 * (x**2 + y**2)
    hz = p.kappa * x
    hx = p.lam * y

    V = np.zeros(x.shape + (2, 2), dtype=float)
    V[..., 0, 0] = common + hz
    V[..., 1, 1] = common - hz
    V[..., 0, 1] = hx
    V[..., 1, 0] = hx
    return V


def adiabatic_energies_2d(x, y, p=LVC2DParameters()):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x, y = np.broadcast_arrays(x, y)

    common = 0.5 * p.omega**2 * (x**2 + y**2)
    rho = np.sqrt((p.kappa*x)**2 + (p.lam*y)**2)

    return np.stack([common-rho, common+rho], axis=-1)


def adiabatic_gradients_2d(R, p=LVC2DParameters(), gap_floor=1e-14):
    """Return gradients[state,coord] away from the exact CI."""
    x, y = np.asarray(R, dtype=float)
    rho = np.sqrt((p.kappa*x)**2 + (p.lam*y)**2)
    if rho < gap_floor:
        raise ValueError("Adiabatic gradient is direction dependent at the exact CI.")

    common = p.omega**2 * np.array([x, y], dtype=float)
    split = np.array(
        [p.kappa**2*x/rho, p.lam**2*y/rho],
        dtype=float,
    )
    return np.stack([common-split, common+split])


def mixing_angle(R, p=LVC2DParameters()):
    x, y = np.asarray(R, dtype=float)
    return np.arctan2(p.lam*y, p.kappa*x)


def analytic_adiabatic_vectors(R, p=LVC2DParameters(), gap_floor=1e-14):
    """Columns are lower and upper real adiabatic eigenvectors."""
    x, y = np.asarray(R, dtype=float)
    rho2 = (p.kappa*x)**2 + (p.lam*y)**2
    if rho2 < gap_floor**2:
        raise ValueError("Adiabatic eigenvectors are undefined at the exact CI.")

    th = mixing_angle((x, y), p)
    s = np.sin(0.5*th)
    c = np.cos(0.5*th)

    return np.array([
        [-s, c],
        [ c, s],
    ], dtype=float)


def vector_nac_2d(R, p=LVC2DParameters(), gap_floor=1e-14):
    """Return d[state_i,state_j,coord] in the analytic half-angle gauge."""
    x, y = np.asarray(R, dtype=float)
    denom = (p.kappa*x)**2 + (p.lam*y)**2
    if denom < gap_floor**2:
        raise ValueError("Derivative coupling is singular at the exact CI.")

    a = 0.5 * p.kappa * p.lam * np.array([-y, x]) / denom

    d = np.zeros((2, 2, 2), dtype=float)
    d[0, 1] = a
    d[1, 0] = -a
    return d


def branching_plane_vectors(p=LVC2DParameters()):
    g = np.array([p.kappa, 0.0])
    h = np.array([0.0, p.lam])
    return g, h


def circle_path(radius=1.0, n=1001):
    phi = np.linspace(0.0, 2.0*np.pi, n)
    return np.column_stack([radius*np.cos(phi), radius*np.sin(phi)])


def berry_line_integral(path, p=LVC2DParameters()):
    """Midpoint line integral of d_01 dot dR along a closed path."""
    path = np.asarray(path, dtype=float)
    total = 0.0
    for a, b in zip(path[:-1], path[1:]):
        mid = 0.5*(a+b)
        dR = b-a
        d01 = vector_nac_2d(mid, p)[0, 1]
        total += np.dot(d01, dR)
    return float(total)


def parallel_transport_real_state(path, state=0, p=LVC2DParameters()):
    """Sign-align numerical eigenvectors successively along a path."""
    path = np.asarray(path, dtype=float)
    transported = []

    ref = None
    for R in path:
        V = diabatic_potential_2d(R[0], R[1], p)
        _, U = np.linalg.eigh(V)
        v = U[:, state].copy()

        if ref is not None and np.dot(ref, v) < 0.0:
            v *= -1.0

        transported.append(v)
        ref = v

    transported = np.asarray(transported)
    return transported, float(np.dot(transported[0], transported[-1]))
