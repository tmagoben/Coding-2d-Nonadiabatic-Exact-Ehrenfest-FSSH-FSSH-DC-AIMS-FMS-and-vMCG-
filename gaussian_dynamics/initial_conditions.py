from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class GaussianWignerEnsemble:
    q: np.ndarray
    p: np.ndarray
    seed: int

    def __post_init__(self):
        q = np.asarray(self.q, dtype=float)
        p = np.asarray(self.p, dtype=float)
        if q.ndim != 2 or p.shape != q.shape:
            raise ValueError("q and p must both have shape (nsample, ndim).")
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "p", p)

    @property
    def nsamples(self):
        return self.q.shape[0]

    @property
    def dimension(self):
        return self.q.shape[1]


def gaussian_wigner_covariances(A):
    """Return coordinate and momentum covariance matrices for a pure Gaussian.

    For
        g(q) ~ exp[-1/2 (q-q0)^T A (q-q0) + i p0^T(q-q0)]
    with real symmetric positive-definite A and hbar=1,

        Cov(q) = 1/2 A^-1,
        Cov(p) = 1/2 A.
    """
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be square.")
    if not np.allclose(A, A.T, atol=1e-12):
        raise ValueError("A must be symmetric.")
    if np.min(np.linalg.eigvalsh(A)) <= 0.0:
        raise ValueError("A must be positive definite.")

    return 0.5 * np.linalg.inv(A), 0.5 * A


def sample_gaussian_wigner(q0, p0, A, nsamples, seed=12345):
    """Sample phase-space points from the Wigner distribution of a frozen Gaussian."""
    q0 = np.asarray(q0, dtype=float)
    p0 = np.asarray(p0, dtype=float)
    if q0.shape != p0.shape or q0.ndim != 1:
        raise ValueError("q0 and p0 must be equal-length vectors.")
    if int(nsamples) <= 0:
        raise ValueError("nsamples must be positive.")

    cov_q, cov_p = gaussian_wigner_covariances(A)
    if cov_q.shape != (len(q0), len(q0)):
        raise ValueError("A dimension is incompatible with q0/p0.")

    rng = np.random.default_rng(int(seed))
    q = rng.multivariate_normal(q0, cov_q, size=int(nsamples))
    p = rng.multivariate_normal(p0, cov_p, size=int(nsamples))

    return GaussianWignerEnsemble(q=q, p=p, seed=int(seed))
