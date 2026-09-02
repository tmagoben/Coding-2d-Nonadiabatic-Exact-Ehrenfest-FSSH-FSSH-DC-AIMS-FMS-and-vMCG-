import numpy as np


def harmonic(x, mass=1.0, omega=1.0, x0=0.0):
    y = np.asarray(x) - x0
    return 0.5 * mass * omega**2 * y**2


def harmonic_gradient(x, mass=1.0, omega=1.0, x0=0.0):
    return mass * omega**2 * (np.asarray(x) - x0)


def harmonic_hessian(x, mass=1.0, omega=1.0, x0=0.0):
    x = np.asarray(x)
    return np.zeros_like(x, dtype=float) + mass * omega**2


def quartic(x, k2=0.5, k4=0.02):
    x = np.asarray(x)
    return 0.5 * k2 * x**2 + k4 * x**4


def quartic_gradient(x, k2=0.5, k4=0.02):
    x = np.asarray(x)
    return k2 * x + 4.0 * k4 * x**3


def quartic_hessian(x, k2=0.5, k4=0.02):
    x = np.asarray(x)
    return k2 + 12.0 * k4 * x**2


def double_well(x, a=0.02, b=2.0):
    x = np.asarray(x)
    return a * (x**2 - b**2) ** 2


def double_well_gradient(x, a=0.02, b=2.0):
    x = np.asarray(x)
    return 4.0 * a * x * (x**2 - b**2)


def avoided_crossing_diabatic(x, k=0.02, shift=2.0, delta=0.0, coupling=0.01, beta=1.0):
    """Simple two-state diabatic avoided-crossing model.

    V11 and V22 are displaced harmonic diabats and V12 is a localized Gaussian
    coupling.
    """
    x = np.asarray(x)
    v11 = 0.5 * k * (x + shift) ** 2
    v22 = delta + 0.5 * k * (x - shift) ** 2
    v12 = coupling * np.exp(-beta * x**2)

    if x.ndim == 0:
        return np.array([[v11, v12], [v12, v22]], dtype=float)

    out = np.zeros((len(x), 2, 2), dtype=float)
    out[:, 0, 0] = v11
    out[:, 1, 1] = v22
    out[:, 0, 1] = v12
    out[:, 1, 0] = v12
    return out


def avoided_crossing_diagonal_gradients(x, k=0.02, shift=2.0):
    """Gradients of V11 and V22 only."""
    x = np.asarray(x)
    return np.array([
        k * (x + shift),
        k * (x - shift),
    ])
