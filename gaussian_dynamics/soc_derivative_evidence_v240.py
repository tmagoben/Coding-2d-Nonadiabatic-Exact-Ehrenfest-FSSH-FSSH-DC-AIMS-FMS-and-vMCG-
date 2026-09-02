"""Transported displaced-geometry derivative evidence for v0.24.0."""

from dataclasses import asdict, dataclass

import numpy as np

from .finite_manifold_transport_v233 import certified_transport_from_overlap_v233
from .manifold_transport_v233 import audit_complete_manifold_transport_v233
from .openmolcas_rassi_snapshot_v240 import ParsedOpenMolcasBundleV240


@dataclass(frozen=True)
class ExternalSOCDerivativePolicyV240:
    geometry_tolerance_bohr: float = 1.0e-12
    hermiticity_tolerance: float = 1.0e-10
    spin_free_convergence_tolerance: float = 1.0e-6
    soc_convergence_tolerance: float = 1.0e-7

    def validate(self):
        for name in self.__dataclass_fields__:
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        return self

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ExternalSOCDerivativeEvidenceV240:
    protocol_fingerprint: str
    bundle_fingerprint: str
    K_spin_free: np.ndarray
    K_soc: np.ndarray
    spin_free_step_changes: np.ndarray
    soc_step_changes: np.ndarray
    maximum_hermiticity_residual: float
    minimum_overlap_singular_value: float
    checks: dict
    passed: bool
    policy: dict

    def as_dict(self):
        return {
            "protocol_fingerprint": self.protocol_fingerprint,
            "bundle_fingerprint": self.bundle_fingerprint,
            "K_spin_free_real": self.K_spin_free.real.tolist(),
            "K_spin_free_imag": self.K_spin_free.imag.tolist(),
            "K_soc_real": self.K_soc.real.tolist(),
            "K_soc_imag": self.K_soc.imag.tolist(),
            "spin_free_step_changes": self.spin_free_step_changes.tolist(),
            "soc_step_changes": self.soc_step_changes.tolist(),
            "maximum_hermiticity_residual": float(
                self.maximum_hermiticity_residual
            ),
            "minimum_overlap_singular_value": float(
                self.minimum_overlap_singular_value
            ),
            "checks": dict(self.checks),
            "passed": bool(self.passed),
            "policy": dict(self.policy),
        }


def _hermiticity_residual_v240(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    return float(np.linalg.norm(matrix - matrix.conj().T, ord="fro"))


def _time_reversal_residual_v240(matrix, J):
    matrix = np.asarray(matrix, dtype=complex)
    return float(np.linalg.norm(matrix - J @ matrix.conj() @ J.conj().T, ord="fro"))


def _spin_free_root_residual_v240(matrix, projectors):
    matrix = np.asarray(matrix, dtype=complex)
    identity = np.eye(matrix.shape[0], dtype=complex)
    residuals = []
    for projector in projectors.values():
        projector = np.asarray(projector, dtype=complex)
        rank = int(round(float(np.trace(projector).real)))
        scalar = np.trace(projector @ matrix) / rank
        residuals.extend(
            (
                np.linalg.norm(projector @ matrix @ projector - scalar * projector, ord="fro"),
                np.linalg.norm(projector @ matrix @ (identity - projector), ord="fro"),
            )
        )
    return float(max(residuals, default=float("inf")))


def audit_external_soc_derivatives_v240(
    bundle,
    *,
    policy=ExternalSOCDerivativePolicyV240(),
):
    """Certify component derivatives after right-to-reference transport.

    The raw overlap is used only to construct its unitary polar transport.  It is
    never substituted directly for transport.
    """

    if type(bundle) is not ParsedOpenMolcasBundleV240:
        raise TypeError("derivative evidence requires a parsed OpenMolcas bundle.")
    policy = policy.validate()
    protocol = bundle.protocol.validate()
    records = bundle.record_map
    if tuple(records) != protocol.expected_record_ids():
        raise ValueError("parsed record inventory is not the protocol inventory.")
    reference = records["reference"]
    reference_geometry = np.asarray(protocol.reference_geometry_bohr, dtype=float)
    geometry_complete = bool(
        np.max(np.abs(reference.geometry_bohr - reference_geometry))
        <= policy.geometry_tolerance_bohr
    )
    convergence_complete = bool(all(all(item.convergence.values()) for item in records.values()))
    ncoord = protocol.coordinate_dimension
    nstate = len(protocol.state_order)
    spin_free_estimates = np.empty((ncoord, 3, nstate, nstate), dtype=complex)
    soc_estimates = np.empty_like(spin_free_estimates)
    geometry_pairs = True
    all_transport_ready = True
    all_manifolds_ready = True
    minimum_singular_value = 1.0
    hermiticity_residuals = []
    time_reversal_residuals = []
    spin_free_root_residuals = []

    for item in records.values():
        hermiticity_residuals.extend(
            (_hermiticity_residual_v240(item.H_spin_free),
             _hermiticity_residual_v240(item.H_soc))
        )
    symmetry = protocol.symmetry_contract()
    model_space = protocol.model_space()
    J = symmetry.time_reversal_matrix
    for item in records.values():
        time_reversal_residuals.extend(
            (
                _time_reversal_residual_v240(item.H_spin_free, J),
                _time_reversal_residual_v240(item.H_soc, J),
            )
        )
        spin_free_root_residuals.append(
            _spin_free_root_residual_v240(item.H_spin_free, symmetry.projectors)
        )
    for coordinate in range(ncoord):
        for step_index, step in enumerate(protocol.displacement_steps_bohr):
            minus = records[f"q{coordinate:02d}_h{step_index}_minus"]
            plus = records[f"q{coordinate:02d}_h{step_index}_plus"]
            expected_minus = reference_geometry.reshape(-1).copy()
            expected_plus = reference_geometry.reshape(-1).copy()
            expected_minus[coordinate] -= step
            expected_plus[coordinate] += step
            geometry_pairs = bool(
                geometry_pairs
                and np.max(np.abs(minus.geometry_bohr.reshape(-1) - expected_minus))
                <= policy.geometry_tolerance_bohr
                and np.max(np.abs(plus.geometry_bohr.reshape(-1) - expected_plus))
                <= policy.geometry_tolerance_bohr
            )
            transported = []
            for record in (minus, plus):
                result = certified_transport_from_overlap_v233(
                    record.reference_overlap
                )
                minimum_singular_value = min(
                    minimum_singular_value, result.minimum_singular_value
                )
                all_transport_ready = bool(
                    all_transport_ready and result.trajectory_ready
                )
                manifold = audit_complete_manifold_transport_v233(
                    record.reference_overlap,
                    model_space,
                    symmetry,
                )
                all_manifolds_ready = bool(all_manifolds_ready and manifold.passed)
                U = result.right_to_left_transport
                transported.append(
                    (
                        U @ record.H_spin_free @ U.conj().T,
                        U @ record.H_soc @ U.conj().T,
                    )
                )
            minus_components, plus_components = transported
            spin_free_estimates[coordinate, step_index] = (
                plus_components[0] - minus_components[0]
            ) / (2.0 * step)
            soc_estimates[coordinate, step_index] = (
                plus_components[1] - minus_components[1]
            ) / (2.0 * step)

    spin_free_changes = np.asarray(
        [
            [
                np.linalg.norm(values[index + 1] - values[index], ord="fro")
                for index in range(2)
            ]
            for values in spin_free_estimates
        ],
        dtype=float,
    )
    soc_changes = np.asarray(
        [
            [
                np.linalg.norm(values[index + 1] - values[index], ord="fro")
                for index in range(2)
            ]
            for values in soc_estimates
        ],
        dtype=float,
    )
    maximum_hermiticity = max(hermiticity_residuals, default=float("inf"))
    derivative_hermiticity = max(
        (
            _hermiticity_residual_v240(matrix)
            for matrix in np.concatenate(
                (spin_free_estimates[:, -1], soc_estimates[:, -1]), axis=0
            )
        ),
        default=float("inf"),
    )
    derivative_symmetry = max(
        (
            *(
                _time_reversal_residual_v240(matrix, J)
                for matrix in spin_free_estimates[:, -1]
            ),
            *(
                _time_reversal_residual_v240(matrix, J)
                for matrix in soc_estimates[:, -1]
            ),
            *(
                _spin_free_root_residual_v240(matrix, symmetry.projectors)
                for matrix in spin_free_estimates[:, -1]
            ),
        ),
        default=float("inf"),
    )
    checks = {
        "exact_record_inventory": tuple(records) == protocol.expected_record_ids(),
        "reference_geometry": geometry_complete,
        "complete_centered_geometry_pairs": geometry_pairs,
        "all_calculations_converged": convergence_complete,
        "component_hermiticity": maximum_hermiticity <= policy.hermiticity_tolerance,
        "derivative_hermiticity": derivative_hermiticity <= policy.hermiticity_tolerance,
        "component_time_reversal_and_spin_degeneracy": bool(
            max(time_reversal_residuals, default=float("inf"))
            <= policy.hermiticity_tolerance
            and max(spin_free_root_residuals, default=float("inf"))
            <= policy.hermiticity_tolerance
        ),
        "derivative_time_reversal_and_spin_degeneracy": (
            derivative_symmetry <= policy.hermiticity_tolerance
        ),
        "finite_manifold_transport": all_transport_ready,
        "complete_manifold_tracking": all_manifolds_ready,
        "spin_free_derivative_convergence": bool(
            np.max(spin_free_changes[:, -1])
            <= policy.spin_free_convergence_tolerance
        ),
        "soc_derivative_convergence": bool(
            np.max(soc_changes[:, -1]) <= policy.soc_convergence_tolerance
        ),
        "finite_spin_free_derivatives": bool(np.all(np.isfinite(spin_free_estimates))),
        "finite_soc_derivatives": bool(np.all(np.isfinite(soc_estimates))),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return ExternalSOCDerivativeEvidenceV240(
        protocol_fingerprint=protocol.fingerprint(),
        bundle_fingerprint=bundle.fingerprint,
        K_spin_free=spin_free_estimates[:, -1].copy(),
        K_soc=soc_estimates[:, -1].copy(),
        spin_free_step_changes=spin_free_changes,
        soc_step_changes=soc_changes,
        maximum_hermiticity_residual=maximum_hermiticity,
        minimum_overlap_singular_value=minimum_singular_value,
        checks=checks,
        passed=bool(all(checks.values())),
        policy=policy.as_dict(),
    )


def require_external_soc_derivatives_v240(*args, **kwargs):
    evidence = audit_external_soc_derivatives_v240(*args, **kwargs)
    if not evidence.passed:
        failed = ", ".join(name for name, value in evidence.checks.items() if not value)
        raise ValueError("external SOC derivative evidence failed: " + failed)
    return evidence
