import numpy as np


def frozen_gaussian(x, q=0.0, p=0.0, alpha=1.0, chirp=0.0):
    """Normalized frozen/generalized Gaussian in atomic units.

    g = (alpha/pi)^(1/4)
        exp[-alpha (x-q)^2/2 + i p(x-q) + i chirp (x-q)^2/2]

    alpha > 0.  chirp is a real quadratic phase.
    """
    if alpha <= 0.0:
        raise ValueError("alpha must be positive.")

    x = np.asarray(x, dtype=float)
    y = x - q
    N = (alpha / np.pi) ** 0.25

    return N * np.exp(
        -0.5 * alpha * y**2
        + 1j * p * y
        + 0.5j * chirp * y**2
    )


def analytic_overlap(qi, pi, qj, pj, alpha):
    """Analytic overlap for equal-width, zero-chirp frozen Gaussians."""
    dq = qi - qj
    dp = pi - pj

    return np.exp(
        -0.25 * alpha * dq**2
        -0.25 * dp**2 / alpha
        + 0.5j * (pi + pj) * dq
    )


def kinetic_on_gaussian(x, q, p, alpha, mass=1.0, chirp=0.0):
    """Apply -1/(2m) d^2/dx^2 analytically to a Gaussian.

    The generalized complex exponent coefficient is z = alpha - i*chirp:
        g ~ exp[-z (x-q)^2/2 + i p(x-q)]
    """
    x = np.asarray(x, dtype=float)
    y = x - q
    g = frozen_gaussian(x, q, p, alpha, chirp)

    z = alpha - 1j * chirp
    f = -z * y + 1j * p
    second_derivative_factor = f**2 - z

    return -(second_derivative_factor * g) / (2.0 * mass)


def gaussian_parameter_derivatives(x, q, p, alpha):
    """Analytic derivatives for the zero-chirp normalized Gaussian."""
    g = frozen_gaussian(x, q, p, alpha)
    y = np.asarray(x) - q

    dq = (alpha * y - 1j * p) * g
    dp = 1j * y * g
    dlogalpha = (0.25 - 0.5 * alpha * y**2) * g

    return dq, dp, dlogalpha


def gaussian_moments(x, psi, dx):
    """Return norm, <x>, variance(x), and <p> using FFT differentiation."""
    psi = np.asarray(psi, dtype=complex)
    nrm = np.vdot(psi, psi).real * dx

    xmean = (np.vdot(psi, x * psi) * dx / nrm).real
    x2mean = (np.vdot(psi, x**2 * psi) * dx / nrm).real
    xvar = x2mean - xmean**2

    k = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    dpsi = np.fft.ifft(1j * k * np.fft.fft(psi))
    pmean = (np.vdot(psi, -1j * dpsi) * dx / nrm).real

    return {
        "norm": nrm,
        "x_mean": xmean,
        "x_variance": xvar,
        "p_mean": pmean,
    }
