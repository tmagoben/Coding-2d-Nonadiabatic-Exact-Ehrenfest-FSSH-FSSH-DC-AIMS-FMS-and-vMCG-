"""Physical-SOC, time-reversal, projector, and force validation for v0.22.0."""

from dataclasses import asdict, dataclass
import numpy as np

from .provider_differential_audit_v214 import (
    ProviderDifferentialAuditSettingsV214,
    audit_provider_differentials_v214,
)
from .gauge_graph import nearest_unitary
from .soc_admission_v221 import (
    audit_soc_symmetry_contract_v221,
    soc_symmetry_contract_from_provider_v221,
)


def _scaled_frobenius_error_v220(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        raise ValueError("SOC validation arrays have incompatible shapes.")
    absolute = float(np.linalg.norm(left - right))
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return absolute, absolute / scale


def time_reversal_residual_v220(matrix, unitary_part):
    matrix = np.asarray(matrix, dtype=complex)
    J = np.asarray(unitary_part, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("time-reversal validation requires a square matrix.")
    if J.shape != matrix.shape:
        raise ValueError("time-reversal matrix has incompatible dimension.")
    transformed = J @ matrix.conj() @ J.conj().T
    return _scaled_frobenius_error_v220(matrix, transformed)[1]


def time_reversal_square_residual_v220(unitary_part, fermionic):
    J = np.asarray(unitary_part, dtype=complex)
    if J.ndim != 2 or J.shape[0] != J.shape[1]:
        raise ValueError("time-reversal unitary part must be square.")
    target = (-1.0 if fermionic else 1.0) * np.eye(J.shape[0], dtype=complex)
    return _scaled_frobenius_error_v220(J @ J.conj(), target)[1]


def transform_time_reversal_matrix_v220(unitary_part, gauge_matrix):
    """Transform J for c'=G^dagger c when Theta=J K."""
    J = np.asarray(unitary_part, dtype=complex)
    G = np.asarray(gauge_matrix, dtype=complex)
    if J.shape != G.shape or G.ndim != 2 or G.shape[0] != G.shape[1]:
        raise ValueError("time-reversal and gauge matrices must be equal-size squares.")
    if not np.all(np.isfinite(J)) or not np.all(np.isfinite(G)):
        raise ValueError("time-reversal and gauge matrices must be finite.")
    if _scaled_frobenius_error_v220(
        G.conj().T @ G, np.eye(G.shape[0], dtype=complex)
    )[1] > 1.0e-12:
        raise ValueError("electronic gauge matrix must be unitary.")
    return G.conj().T @ J @ G.conj()


def transform_projector_v220(projector, gauge_matrix):
    projector = np.asarray(projector, dtype=complex)
    G = np.asarray(gauge_matrix, dtype=complex)
    if projector.shape != G.shape:
        raise ValueError("projector and gauge matrix dimensions differ.")
    if not np.all(np.isfinite(projector)) or not np.all(np.isfinite(G)):
        raise ValueError("projector and gauge matrix must be finite.")
    if _scaled_frobenius_error_v220(
        G.conj().T @ G, np.eye(G.shape[0], dtype=complex)
    )[1] > 1.0e-12:
        raise ValueError("electronic gauge matrix must be unitary.")
    return G.conj().T @ projector @ G


def projector_population_v220(electronic_vector, projector):
    vector = np.asarray(electronic_vector, dtype=complex)
    projector = np.asarray(projector, dtype=complex)
    if vector.ndim != 1 or projector.shape != (len(vector), len(vector)):
        raise ValueError("electronic vector and projector dimensions differ.")
    if not np.all(np.isfinite(vector)) or not np.all(np.isfinite(projector)):
        raise ValueError("electronic vector and projector must be finite.")
    projector_residual = max(
        _scaled_frobenius_error_v220(projector, projector.conj().T)[1],
        _scaled_frobenius_error_v220(projector @ projector, projector)[1],
    )
    if projector_residual > 1.0e-10:
        raise ValueError("population operator must be a Hermitian projector.")
    norm = float(np.real(np.vdot(vector, vector)))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("electronic vector must have finite nonzero norm.")
    population = float(np.real(np.vdot(vector, projector @ vector)) / norm)
    if population < -1.0e-10 or population > 1.0 + 1.0e-10:
        raise ValueError("projector population lies outside [0,1].")
    return float(np.clip(population, 0.0, 1.0))


def _projector_residuals_v220(projectors):
    projectors = tuple(np.asarray(item, dtype=complex) for item in projectors.values())
    if not projectors:
        raise ValueError("at least one physical projector is required.")
    nstate = projectors[0].shape[0]
    residuals = []
    for index, projector in enumerate(projectors):
        if projector.shape != (nstate, nstate):
            raise ValueError("physical projectors have incompatible dimensions.")
        residuals.append(_scaled_frobenius_error_v220(projector, projector.conj().T)[1])
        residuals.append(_scaled_frobenius_error_v220(projector @ projector, projector)[1])
        for other in projectors[index + 1 :]:
            residuals.append(float(np.linalg.norm(projector @ other)))
    residuals.append(
        _scaled_frobenius_error_v220(
            sum(projectors, np.zeros_like(projectors[0])), np.eye(nstate)
        )[1]
    )
    return max(residuals, default=0.0)


@dataclass(frozen=True)
class PhysicalSOCAuditSettingsV220:
    composition_tolerance: float = 1.0e-13
    time_reversal_tolerance: float = 1.0e-12
    projector_tolerance: float = 1.0e-12
    force_tolerance: float = 2.0e-10
    force_difference_step: float = 1.0e-5
    component_difference_steps: tuple[float, ...] = (1.0e-4, 5.0e-5, 2.5e-5)
    component_derivative_tolerance: float = 2.0e-9
    differential_settings: ProviderDifferentialAuditSettingsV214 = (
        ProviderDifferentialAuditSettingsV214()
    )

    def validate(self):
        for name, value in (
            ("composition_tolerance", self.composition_tolerance),
            ("time_reversal_tolerance", self.time_reversal_tolerance),
            ("projector_tolerance", self.projector_tolerance),
            ("force_tolerance", self.force_tolerance),
            ("force_difference_step", self.force_difference_step),
            ("component_derivative_tolerance", self.component_derivative_tolerance),
        ):
            if not np.isfinite(value) or float(value) <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        self.differential_settings.validate()
        steps = np.asarray(self.component_difference_steps, dtype=float)
        if (
            steps.ndim != 1
            or len(steps) < 2
            or not np.all(np.isfinite(steps))
            or np.any(steps <= 0.0)
        ):
            raise ValueError(
                "component_difference_steps must contain at least two positive values."
            )
        return self


@dataclass(frozen=True)
class ComponentDerivativeRowV221:
    coordinate: int
    step: float
    spin_free_scaled_error: float
    soc_scaled_error: float

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class PhysicalSOCAuditV220:
    q: np.ndarray
    H_composition_error: float
    K_composition_error: float
    maximum_time_reversal_residual: float
    time_reversal_square_residual: float
    projector_residual: float
    soc_force_error: float
    soc_force_analytic: float
    soc_force_finite_difference: float
    maximum_spin_free_component_derivative_error: float
    maximum_soc_component_derivative_error: float
    component_derivative_rows: tuple[ComponentDerivativeRowV221, ...]
    symmetry_report: object
    differential_report: object
    provenance_fingerprint: str
    checks: dict
    passed: bool
    thresholds: dict

    def as_dict(self):
        return {
            "q": np.asarray(self.q, dtype=float).tolist(),
            "H_composition_error": self.H_composition_error,
            "K_composition_error": self.K_composition_error,
            "maximum_time_reversal_residual": self.maximum_time_reversal_residual,
            "time_reversal_square_residual": self.time_reversal_square_residual,
            "projector_residual": self.projector_residual,
            "soc_force_error": self.soc_force_error,
            "soc_force_analytic": self.soc_force_analytic,
            "soc_force_finite_difference": self.soc_force_finite_difference,
            "maximum_spin_free_component_derivative_error": (
                self.maximum_spin_free_component_derivative_error
            ),
            "maximum_soc_component_derivative_error": (
                self.maximum_soc_component_derivative_error
            ),
            "component_derivative_rows": [
                row.as_dict() for row in self.component_derivative_rows
            ],
            "symmetry_report": self.symmetry_report.as_dict(),
            "differential_report": self.differential_report.as_dict(),
            "provenance_fingerprint": self.provenance_fingerprint,
            "checks": dict(self.checks),
            "passed": bool(self.passed),
            "thresholds": dict(self.thresholds),
        }


def audit_physical_soc_provider_v220(
    provider,
    q,
    *,
    fermionic=None,
    settings=PhysicalSOCAuditSettingsV220(),
):
    """Certify decomposition, derivatives, time reversal, projectors, and SOC force."""
    settings = settings.validate()
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or len(q) < 1 or not np.all(np.isfinite(q)):
        raise ValueError("physical SOC audit requires a finite coordinate vector.")
    provenance = provider.provenance.validate()
    components = provider.components(q).validate()
    snapshot = provider.evaluate_snapshot(q).validate()
    fingerprint = provenance.fingerprint()

    _, H_error = _scaled_frobenius_error_v220(snapshot.point.H, components.H)
    _, K_error = _scaled_frobenius_error_v220(
        snapshot.point.hamiltonian_derivative_operator_q, components.K
    )
    symmetry_contract = soc_symmetry_contract_from_provider_v221(provider)
    symmetry_report = audit_soc_symmetry_contract_v221(
        provenance.model_space,
        symmetry_contract,
        provenance=provenance,
        fermionic=fermionic,
        tolerance=settings.time_reversal_tolerance,
    )
    J = symmetry_contract.time_reversal_matrix
    time_reversal_residuals = [
        time_reversal_residual_v220(matrix, J)
        for matrix in (
            components.H_spin_free,
            components.H_soc,
            components.H,
            *components.K_spin_free,
            *components.K_soc,
            *components.K,
        )
    ]
    maximum_time_reversal = max(time_reversal_residuals, default=0.0)
    square_residual = symmetry_report.time_reversal_square_residual
    projector_residual = symmetry_report.projector_residual

    component_rows = []
    for coordinate in range(len(q)):
        for step in settings.component_difference_steps:
            displacement = np.zeros_like(q)
            displacement[coordinate] = float(step)
            plus_snapshot = provider.evaluate_snapshot(q + displacement).validate()
            minus_snapshot = provider.evaluate_snapshot(q - displacement).validate()
            plus_components = provider.components(q + displacement).validate()
            minus_components = provider.components(q - displacement).validate()
            overlap_plus = np.asarray(
                provider.snapshot_overlap(snapshot, plus_snapshot), dtype=complex
            )
            overlap_minus = np.asarray(
                provider.snapshot_overlap(snapshot, minus_snapshot), dtype=complex
            )
            if overlap_plus.shape != components.H.shape or overlap_minus.shape != components.H.shape:
                raise ValueError("component differential overlap has incompatible shape.")
            if not np.all(np.isfinite(overlap_plus)) or not np.all(np.isfinite(overlap_minus)):
                raise ValueError("component differential overlap contains non-finite data.")
            unitary_plus = nearest_unitary(overlap_plus)
            unitary_minus = nearest_unitary(overlap_minus)

            def centered(name):
                plus_matrix = getattr(plus_components, name)
                minus_matrix = getattr(minus_components, name)
                plus_center = unitary_plus @ plus_matrix @ unitary_plus.conj().T
                minus_center = unitary_minus @ minus_matrix @ unitary_minus.conj().T
                return (plus_center - minus_center) / (2.0 * float(step))

            finite_spin_free = centered("H_spin_free")
            finite_soc = centered("H_soc")
            _, spin_free_error = _scaled_frobenius_error_v220(
                finite_spin_free, components.K_spin_free[coordinate]
            )
            _, soc_error = _scaled_frobenius_error_v220(
                finite_soc, components.K_soc[coordinate]
            )
            component_rows.append(
                ComponentDerivativeRowV221(
                    coordinate=int(coordinate),
                    step=float(step),
                    spin_free_scaled_error=float(spin_free_error),
                    soc_scaled_error=float(soc_error),
                )
            )
    maximum_spin_free_component = max(
        (row.spin_free_scaled_error for row in component_rows), default=0.0
    )
    maximum_soc_component = max(
        (row.soc_scaled_error for row in component_rows), default=0.0
    )

    indices = np.arange(1, components.H.shape[0] + 1, dtype=float)
    vector = (1.0 + 0.13 * indices) + 1j * (0.19 - 0.07 * indices)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    step = float(settings.force_difference_step)
    force_rows = []
    for coordinate in range(len(q)):
        displacement = np.zeros_like(q)
        displacement[coordinate] = step
        plus = provider.components(q + displacement).H_soc
        minus = provider.components(q - displacement).H_soc
        analytic = -float(
            np.real(np.trace(density @ components.K_soc[coordinate]))
        )
        plus_energy = float(np.real(np.trace(density @ plus)))
        minus_energy = float(np.real(np.trace(density @ minus)))
        finite = -(plus_energy - minus_energy) / (2.0 * step)
        force_rows.append((analytic, finite, abs(analytic - finite)))
    analytic_force, finite_force, _ = force_rows[0]
    force_error = max((row[2] for row in force_rows), default=0.0)

    differential = audit_provider_differentials_v214(
        provider,
        q,
        provenance,
        settings=settings.differential_settings,
    )
    emitted_fingerprint = snapshot.point.metadata.get("v213_provenance_fingerprint")
    soc_expected = bool(
        provenance.parameters.get("soc_signal_expected", provenance.soc_enabled)
    )
    soc_signal = max(
        float(np.linalg.norm(components.H_soc)),
        float(np.linalg.norm(components.K_soc)),
    )
    checks = {
        "H_decomposition": H_error <= settings.composition_tolerance,
        "K_decomposition": K_error <= settings.composition_tolerance,
        "time_reversal": maximum_time_reversal <= settings.time_reversal_tolerance,
        "SOC_symmetry_admission": symmetry_report.passed,
        "time_reversal_unitarity": symmetry_report.checks[
            "time_reversal_unitarity"
        ],
        "time_reversal_square": square_residual <= settings.time_reversal_tolerance,
        "physical_projectors": projector_residual <= settings.projector_tolerance,
        "spin_free_component_derivatives": maximum_spin_free_component
        <= settings.component_derivative_tolerance,
        "SOC_component_derivatives": maximum_soc_component
        <= settings.component_derivative_tolerance,
        "SOC_force_derivative": force_error <= settings.force_tolerance,
        "cross_geometry_differentials": differential.passed,
        "provider_provenance": emitted_fingerprint == fingerprint,
        "physical_SOC_signal": (soc_signal > 0.0) if soc_expected else (soc_signal == 0.0),
    }
    return PhysicalSOCAuditV220(
        q=q.copy(),
        H_composition_error=float(H_error),
        K_composition_error=float(K_error),
        maximum_time_reversal_residual=float(maximum_time_reversal),
        time_reversal_square_residual=float(square_residual),
        projector_residual=float(projector_residual),
        soc_force_error=float(force_error),
        soc_force_analytic=float(analytic_force),
        soc_force_finite_difference=float(finite_force),
        maximum_spin_free_component_derivative_error=float(
            maximum_spin_free_component
        ),
        maximum_soc_component_derivative_error=float(maximum_soc_component),
        component_derivative_rows=tuple(component_rows),
        symmetry_report=symmetry_report,
        differential_report=differential,
        provenance_fingerprint=fingerprint,
        checks=checks,
        passed=bool(all(checks.values())),
        thresholds=asdict(settings),
    )


def require_physical_soc_contract_v220(*args, **kwargs):
    report = audit_physical_soc_provider_v220(*args, **kwargs)
    if not report.passed:
        failed = ", ".join(name for name, passed in report.checks.items() if not passed)
        raise ValueError(f"physical SOC contract failed: {failed}.")
    return report


@dataclass(frozen=True)
class KramersAuditV220:
    geometries: int
    maximum_H_time_reversal_residual: float
    maximum_K_time_reversal_residual: float
    maximum_pair_splitting: float
    time_reversal_square_residual: float
    passed: bool
    tolerance: float

    def as_dict(self):
        return asdict(self)


def audit_kramers_degeneracy_v220(provider, geometries, *, tolerance=1.0e-11):
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("Kramers tolerance must be finite and positive.")
    geometries = tuple(np.asarray(q, dtype=float) for q in geometries)
    if not geometries:
        raise ValueError("Kramers audit requires at least one geometry.")
    contract = soc_symmetry_contract_from_provider_v221(provider)
    symmetry = audit_soc_symmetry_contract_v221(
        provider.provenance.model_space,
        contract,
        provenance=provider.provenance,
        fermionic=True,
        tolerance=tolerance,
    )
    J = contract.time_reversal_matrix
    H_residuals = []
    K_residuals = []
    splittings = []
    for q in geometries:
        point = provider.evaluate_snapshot(q).point
        H_residuals.append(time_reversal_residual_v220(point.H, J))
        for derivative in point.hamiltonian_derivative_operator_q:
            K_residuals.append(time_reversal_residual_v220(derivative, J))
        energies = np.linalg.eigvalsh(point.H)
        if len(energies) % 2:
            raise ValueError("Kramers model dimension must be even.")
        splittings.extend(abs(energies[1::2] - energies[0::2]).tolist())
    maximum_H = max(H_residuals, default=0.0)
    maximum_K = max(K_residuals, default=0.0)
    maximum_splitting = max(splittings, default=0.0)
    square = time_reversal_square_residual_v220(J, fermionic=True)
    passed = (
        symmetry.passed
        and max(maximum_H, maximum_K, maximum_splitting, square) <= tolerance
    )
    return KramersAuditV220(
        geometries=len(geometries),
        maximum_H_time_reversal_residual=float(maximum_H),
        maximum_K_time_reversal_residual=float(maximum_K),
        maximum_pair_splitting=float(maximum_splitting),
        time_reversal_square_residual=float(square),
        passed=bool(passed),
        tolerance=tolerance,
    )
