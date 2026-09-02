import numpy as np


def uniform_grid(xmin=-12.0, xmax=12.0, n=2048):
    """Periodic FFT grid with the endpoint excluded."""
    x = np.linspace(xmin, xmax, n, endpoint=False)
    dx = (xmax - xmin) / n
    return x, dx


def inner_product(a, b, dx):
    """Discrete approximation to integral a*(x)b(x) dx."""
    return np.vdot(a, b) * dx


def norm(psi, dx):
    return float(np.sqrt(max(inner_product(psi, psi, dx).real, 0.0)))


def normalize(psi, dx):
    nrm = norm(psi, dx)
    if nrm == 0.0:
        raise ValueError("Cannot normalize a zero wavefunction.")
    return np.asarray(psi, dtype=complex) / nrm
