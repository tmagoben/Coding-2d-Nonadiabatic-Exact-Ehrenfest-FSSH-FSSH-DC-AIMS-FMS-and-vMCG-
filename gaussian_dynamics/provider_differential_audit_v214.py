"""Cross-geometry certification of the frozen H/K/D electronic contract.

Pointwise Hermiticity is necessary but not sufficient.  A provider can return
individually valid H, K, and D matrices that do not describe one differentiable
electronic model.  v0.21.4 therefore compares analytic K and D against centered
finite differences built from cross-geometry electronic overlaps.
"""

from dataclasses import asdict, dataclass
import numpy as np

from .electronic_contract_v213 import validate_electronic_contract_v213
from .gauge_graph import nearest_unitary
from .matrix_invariants_v213 import (
    antihermiticity_residual_v213,
    hermiticity_residual_v213,
    isometry_residual_v213,
    symmetry_residual_v213,
)


def _scaled_frobenius_error(left, right, scale_floor=1.0):
    a = np.asarray(left, dtype=complex)
    b = np.asarray(right, dtype=complex)
    if a.shape != b.shape:
        raise ValueError("differential-audit matrices have incompatible shapes.")
    absolute = float(np.linalg.norm(a - b, ord="fro"))
    scale = max(
        float(np.linalg.norm(a, ord="fro")),
        float(np.linalg.norm(b, ord="fro")),
        float(scale_floor),
    )
    return absolute, absolute / scale


@dataclass(frozen=True)
class ProviderDifferentialAuditSettingsV214:
    coordinate_steps: tuple[float, ...] | None = None
    default_step: float = 1.0e-4
    structural_tolerance: float = 1.0e-12
    hamiltonian_derivative_tolerance: float = 2.0e-9
    connection_tolerance: float = 2.0e-9
    overlap_isometry_tolerance: float = 1.0e-8

    def steps_for_dimension(self, dimension):
        dimension = int(dimension)
        if dimension < 1:
            raise ValueError("the differential audit needs at least one coordinate.")
        if self.coordinate_steps is None:
            steps = np.full(dimension, float(self.default_step), dtype=float)
        else:
            steps = np.asarray(self.coordinate_steps, dtype=float)
            if steps.shape != (dimension,):
                raise ValueError(
                    "coordinate_steps must contain one value per nuclear coordinate."
                )
        if not np.all(np.isfinite(steps)) or np.any(steps <= 0.0):
            raise ValueError("all finite-difference steps must be finite and positive.")
        return steps

    def validate(self):
        for name, value in (
            ("default_step", self.default_step),
            ("structural_tolerance", self.structural_tolerance),
            (
                "hamiltonian_derivative_tolerance",
                self.hamiltonian_derivative_tolerance,
            ),
            ("connection_tolerance", self.connection_tolerance),
            ("overlap_isometry_tolerance", self.overlap_isometry_tolerance),
        ):
            if not np.isfinite(value) or float(value) <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.coordinate_steps is not None:
            raw = np.asarray(self.coordinate_steps, dtype=float)
            if raw.ndim != 1 or not np.all(np.isfinite(raw)) or np.any(raw <= 0.0):
                raise ValueError(
                    "coordinate_steps must be a finite positive one-dimensional vector."
                )
        return self


@dataclass(frozen=True)
class CoordinateDifferentialAuditV214:
    coordinate: int
    step: float
    hamiltonian_derivative_absolute_error: float
    hamiltonian_derivative_scaled_error: float
    connection_absolute_error: float
    connection_scaled_error: float
    plus_overlap_isometry_residual: float
    minus_overlap_isometry_residual: float
    finite_difference_K_hermiticity_residual: float
    finite_difference_D_antihermiticity_residual: float

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ProviderDifferentialAuditV214:
    q: np.ndarray
    rows: tuple[CoordinateDifferentialAuditV214, ...]
    maximum_structural_residual: float
    maximum_hamiltonian_derivative_scaled_error: float
    maximum_connection_scaled_error: float
    maximum_overlap_isometry_residual: float
    provenance_fingerprint: str
    passed: bool
    checks: dict
    thresholds: dict

    def as_dict(self):
        return {
            "q": np.asarray(self.q, dtype=float).tolist(),
            "rows": [row.as_dict() for row in self.rows],
            "maximum_structural_residual": self.maximum_structural_residual,
            "maximum_hamiltonian_derivative_scaled_error": (
                self.maximum_hamiltonian_derivative_scaled_error
            ),
            "maximum_connection_scaled_error": self.maximum_connection_scaled_error,
            "maximum_overlap_isometry_residual": (
                self.maximum_overlap_isometry_residual
            ),
            "provenance_fingerprint": self.provenance_fingerprint,
            "passed": bool(self.passed),
            "checks": dict(self.checks),
            "thresholds": dict(self.thresholds),
        }


def _validated_overlap(provider, left, right, nstate):
    overlap = np.asarray(provider.snapshot_overlap(left, right), dtype=complex)
    if overlap.shape != (nstate, nstate):
        raise ValueError("cross-geometry overlap has incompatible shape.")
    if not np.all(np.isfinite(overlap)):
        raise ValueError("cross-geometry overlap contains non-finite data.")
    return overlap


def _maximum_point_structural_residual(point):
    residuals = point.structural_residuals_v213()
    return max(
        residuals["H_hermiticity"],
        *residuals["dH_hermiticity"],
        *residuals["connection_antihermiticity"],
        residuals["mass_symmetry"],
    )


def audit_provider_differentials_v214(
    provider,
    q,
    provenance,
    *,
    settings=ProviderDifferentialAuditSettingsV214(),
):
    """Audit K and D against overlap-transported centered differences.

    With ``O_0+ = <Phi(q)|Phi(q+h)>``, the plus-point Hamiltonian expressed in
    the center frame is ``U_0+ H(q+h) U_0+^dagger``, where U is the nearest
    unitary polar factor.  Its centered derivative must equal the physical
    derivative operator K, not the naive derivative of a moving-frame H.

    The centered derivative of the *raw* overlap is compared with D.  This is
    gauge covariant and independently detects a provider that supplies a
    plausible H/K pair but an inconsistent connection.
    """
    settings = settings.validate()
    provenance = provenance.validate()
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or not np.all(np.isfinite(q)):
        raise ValueError("audit geometry must be a finite coordinate vector.")
    steps = settings.steps_for_dimension(len(q))
    fingerprint = provenance.fingerprint()

    center = provider.evaluate_snapshot(q).validate()
    validate_electronic_contract_v213(
        center.point, provenance, tolerance=settings.structural_tolerance
    )
    if not np.array_equal(center.point.q, q):
        raise ValueError("provider returned a geometry different from the request.")
    nstate = center.point.nstate
    if center.state_vectors is None and center.wavefunction_snapshot is None:
        raise ValueError(
            "differential certification requires a cross-geometry overlap path."
        )

    rows = []
    structural = [_maximum_point_structural_residual(center.point)]
    fingerprints = {
        center.point.metadata.get("v213_provenance_fingerprint")
    }
    for coordinate, step in enumerate(steps):
        displacement = np.zeros_like(q)
        displacement[coordinate] = step
        plus = provider.evaluate_snapshot(q + displacement).validate()
        minus = provider.evaluate_snapshot(q - displacement).validate()
        for requested, snapshot in ((q + displacement, plus), (q - displacement, minus)):
            validate_electronic_contract_v213(
                snapshot.point,
                provenance,
                tolerance=settings.structural_tolerance,
            )
            if not np.array_equal(snapshot.point.q, requested):
                raise ValueError("provider returned a displaced geometry incorrectly.")
            structural.append(_maximum_point_structural_residual(snapshot.point))
            fingerprints.add(
                snapshot.point.metadata.get("v213_provenance_fingerprint")
            )

        overlap_plus = _validated_overlap(provider, center, plus, nstate)
        overlap_minus = _validated_overlap(provider, center, minus, nstate)
        unitary_plus = nearest_unitary(overlap_plus)
        unitary_minus = nearest_unitary(overlap_minus)

        H_plus_center = unitary_plus @ plus.point.H @ unitary_plus.conj().T
        H_minus_center = unitary_minus @ minus.point.H @ unitary_minus.conj().T
        K_finite_difference = (H_plus_center - H_minus_center) / (2.0 * step)
        D_finite_difference = (overlap_plus - overlap_minus) / (2.0 * step)

        K_absolute, K_scaled = _scaled_frobenius_error(
            K_finite_difference,
            center.point.hamiltonian_derivative_operator_q[coordinate],
        )
        D_absolute, D_scaled = _scaled_frobenius_error(
            D_finite_difference,
            center.point.connection_q[coordinate],
        )
        rows.append(
            CoordinateDifferentialAuditV214(
                coordinate=int(coordinate),
                step=float(step),
                hamiltonian_derivative_absolute_error=K_absolute,
                hamiltonian_derivative_scaled_error=K_scaled,
                connection_absolute_error=D_absolute,
                connection_scaled_error=D_scaled,
                plus_overlap_isometry_residual=isometry_residual_v213(overlap_plus),
                minus_overlap_isometry_residual=isometry_residual_v213(overlap_minus),
                finite_difference_K_hermiticity_residual=(
                    hermiticity_residual_v213(K_finite_difference)
                ),
                finite_difference_D_antihermiticity_residual=(
                    antihermiticity_residual_v213(D_finite_difference)
                ),
            )
        )

    maximum_structural = max(structural, default=0.0)
    maximum_K = max(
        (row.hamiltonian_derivative_scaled_error for row in rows), default=0.0
    )
    maximum_D = max((row.connection_scaled_error for row in rows), default=0.0)
    maximum_overlap = max(
        (
            max(
                row.plus_overlap_isometry_residual,
                row.minus_overlap_isometry_residual,
            )
            for row in rows
        ),
        default=0.0,
    )
    fingerprint_consistent = fingerprints == {fingerprint}
    checks = {
        "structural_invariants": maximum_structural
        <= settings.structural_tolerance,
        "physical_H_derivatives": maximum_K
        <= settings.hamiltonian_derivative_tolerance,
        "derivative_connections": maximum_D <= settings.connection_tolerance,
        "cross_geometry_overlap_isometry": maximum_overlap
        <= settings.overlap_isometry_tolerance,
        "provenance_fingerprint_consistency": fingerprint_consistent,
    }
    return ProviderDifferentialAuditV214(
        q=q.copy(),
        rows=tuple(rows),
        maximum_structural_residual=float(maximum_structural),
        maximum_hamiltonian_derivative_scaled_error=float(maximum_K),
        maximum_connection_scaled_error=float(maximum_D),
        maximum_overlap_isometry_residual=float(maximum_overlap),
        provenance_fingerprint=fingerprint,
        passed=bool(all(checks.values())),
        checks=checks,
        thresholds=asdict(settings),
    )


def require_provider_differential_contract_v214(*args, **kwargs):
    report = audit_provider_differentials_v214(*args, **kwargs)
    if not report.passed:
        failed = ", ".join(name for name, value in report.checks.items() if not value)
        raise ValueError(f"electronic provider differential contract failed: {failed}.")
    return report
