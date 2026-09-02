"""Strict structural residuals for the v0.21.3 SOC-contract freeze.

The historical code often passed ``atol`` to ``numpy.allclose`` without overriding
NumPy's comparatively loose default relative tolerance.  Structural identities such
as Hermiticity and unitarity need an explicit, inspectable residual instead.
"""

from dataclasses import dataclass, asdict
import numpy as np


@dataclass(frozen=True)
class MatrixInvariantTolerancesV213:
    hermiticity: float = 1.0e-12
    antihermiticity: float = 1.0e-12
    isometry: float = 1.0e-10
    mass_symmetry: float = 1.0e-12
    positive_eigenvalue_floor: float = 0.0

    def validate(self):
        values = asdict(self)
        if any(float(value) < 0.0 for value in values.values()):
            raise ValueError("matrix-invariant tolerances cannot be negative.")
        return self


def scaled_matrix_residual_v213(left, right, *, scale_floor=1.0):
    """Return ``||left-right||_F / max(||left||_F, ||right||_F, floor)``."""
    a = np.asarray(left, dtype=complex)
    b = np.asarray(right, dtype=complex)
    if a.shape != b.shape:
        raise ValueError("matrix residual operands must have identical shapes.")
    numerator = np.linalg.norm(a - b, ord="fro")
    denominator = max(
        np.linalg.norm(a, ord="fro"),
        np.linalg.norm(b, ord="fro"),
        float(scale_floor),
    )
    return float(numerator / denominator)


def hermiticity_residual_v213(matrix):
    a = np.asarray(matrix, dtype=complex)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("Hermiticity requires a square matrix.")
    return scaled_matrix_residual_v213(a, a.conj().T)


def antihermiticity_residual_v213(matrix):
    a = np.asarray(matrix, dtype=complex)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("anti-Hermiticity requires a square matrix.")
    return scaled_matrix_residual_v213(a, -a.conj().T)


def isometry_residual_v213(matrix):
    a = np.asarray(matrix, dtype=complex)
    if a.ndim != 2:
        raise ValueError("an isometry must be a matrix.")
    gram = a.conj().T @ a
    return scaled_matrix_residual_v213(gram, np.eye(a.shape[1], dtype=complex))


def symmetry_residual_v213(matrix):
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("symmetry requires a square matrix.")
    numerator = np.linalg.norm(a - a.T, ord="fro")
    denominator = max(np.linalg.norm(a, ord="fro"), 1.0)
    return float(numerator / denominator)


def require_residual_v213(name, residual, tolerance):
    residual = float(residual)
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError(f"{name} tolerance must be finite and nonnegative.")
    if not np.isfinite(residual):
        raise ValueError(f"{name} residual is non-finite.")
    if residual > tolerance:
        raise ValueError(
            f"{name} residual {residual:.6e} exceeds tolerance {tolerance:.6e}."
        )
    return residual
