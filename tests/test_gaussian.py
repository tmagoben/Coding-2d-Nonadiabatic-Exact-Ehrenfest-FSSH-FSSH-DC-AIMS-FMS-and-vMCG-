import numpy as np

from gaussian_dynamics import (
    uniform_grid,
    frozen_gaussian,
    analytic_overlap,
    inner_product,
    gaussian_moments,
)


def test_normalized_gaussian_and_moments():
    x, dx = uniform_grid(-12.0, 12.0, 4096)
    alpha = 1.4
    q = -0.7
    p = 0.9

    g = frozen_gaussian(x, q, p, alpha)
    m = gaussian_moments(x, g, dx)

    assert abs(m["norm"] - 1.0) < 1e-12
    assert abs(m["x_mean"] - q) < 1e-11
    assert abs(m["x_variance"] - 1.0 / (2.0 * alpha)) < 1e-10
    assert abs(m["p_mean"] - p) < 1e-10


def test_analytic_overlap_matches_quadrature():
    x, dx = uniform_grid(-15.0, 15.0, 4096)

    qi, pi = -1.1, 0.8
    qj, pj = 0.6, -0.2
    alpha = 1.3

    gi = frozen_gaussian(x, qi, pi, alpha)
    gj = frozen_gaussian(x, qj, pj, alpha)

    numeric = inner_product(gi, gj, dx)
    analytic = analytic_overlap(qi, pi, qj, pj, alpha)

    assert abs(numeric - analytic) < 1e-11
