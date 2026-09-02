"""Complete-multiplet and Kramers-safe overlap tracking for v0.23.3."""

from dataclasses import asdict, dataclass

import numpy as np

from .soc_admission_v221 import (
    SOCSymmetryContractV221,
    audit_soc_symmetry_contract_v221,
)


@dataclass(frozen=True)
class ManifoldTransportPolicyV233:
    minimum_assigned_singular_value: float = 0.5
    maximum_competing_leakage: float = 0.2
    minimum_assignment_margin: float = 0.3
    time_reversal_tolerance: float = 1.0e-10
    projector_tolerance: float = 1.0e-10

    def validate(self):
        for name in (
            "minimum_assigned_singular_value",
            "maximum_competing_leakage",
            "minimum_assignment_margin",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0,1].")
        for name in ("time_reversal_tolerance", "projector_tolerance"):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        return self

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ManifoldBlockTransportV233:
    name: str
    dimension: int
    singular_values: np.ndarray
    minimum_singular_value: float
    maximum_competing_leakage: float
    assignment_margin: float

    def as_dict(self):
        return {
            "name": self.name,
            "dimension": int(self.dimension),
            "singular_values": np.asarray(self.singular_values, dtype=float).tolist(),
            "minimum_singular_value": float(self.minimum_singular_value),
            "maximum_competing_leakage": float(
                self.maximum_competing_leakage
            ),
            "assignment_margin": float(self.assignment_margin),
        }


@dataclass(frozen=True)
class CompleteManifoldTransportAuditV233:
    electron_parity: str
    manifold_blocks: tuple
    time_reversal_covariance_residual: float
    checks: dict
    passed: bool
    policy: dict

    def as_dict(self):
        return {
            "electron_parity": self.electron_parity,
            "manifold_blocks": [item.as_dict() for item in self.manifold_blocks],
            "time_reversal_covariance_residual": float(
                self.time_reversal_covariance_residual
            ),
            "checks": dict(self.checks),
            "passed": bool(self.passed),
            "policy": dict(self.policy),
        }


def _projector_basis_v233(projector, *, tolerance):
    projector = np.asarray(projector, dtype=complex)
    hermitian = 0.5 * (projector + projector.conj().T)
    values, vectors = np.linalg.eigh(hermitian)
    selected = values > 0.5
    if not np.any(selected):
        raise ValueError("physical manifold projector has zero rank.")
    if np.max(np.abs(values - np.rint(values))) > tolerance:
        raise ValueError("physical manifold projector is not idempotent.")
    return vectors[:, selected]


def audit_complete_manifold_transport_v233(
    overlap_left_right,
    left_model_space,
    left_symmetry_contract,
    *,
    right_model_space=None,
    right_symmetry_contract=None,
    left_provenance=None,
    right_provenance=None,
    policy=ManifoldTransportPolicyV233(),
):
    """Audit complete physical manifolds in independently gauged endpoint frames."""
    policy = policy.validate()
    left_model_space = left_model_space.validate()
    if right_model_space is None:
        right_model_space = left_model_space
    right_model_space = right_model_space.validate()
    if right_symmetry_contract is None:
        right_symmetry_contract = left_symmetry_contract
    if not isinstance(left_symmetry_contract, SOCSymmetryContractV221) or not isinstance(
        right_symmetry_contract, SOCSymmetryContractV221
    ):
        raise TypeError("manifold transport requires two SOC symmetry contracts.")
    if left_model_space.nstate != right_model_space.nstate:
        raise ValueError("endpoint electronic dimensions differ.")
    if left_symmetry_contract.electron_parity != right_symmetry_contract.electron_parity:
        raise ValueError("endpoint electron-parity sectors differ.")
    if set(left_symmetry_contract.projectors) != set(
        right_symmetry_contract.projectors
    ):
        raise ValueError("endpoint physical-manifold names differ.")

    left_audit = audit_soc_symmetry_contract_v221(
        left_model_space,
        left_symmetry_contract,
        provenance=left_provenance,
        tolerance=policy.projector_tolerance,
    )
    right_audit = audit_soc_symmetry_contract_v221(
        right_model_space,
        right_symmetry_contract,
        provenance=right_provenance,
        tolerance=policy.projector_tolerance,
    )
    overlap = np.asarray(overlap_left_right, dtype=complex)
    nstate = left_model_space.nstate
    if overlap.shape != (nstate, nstate) or not np.all(np.isfinite(overlap)):
        raise ValueError("manifold overlap has incompatible or non-finite data.")

    names = tuple(sorted(left_symmetry_contract.projectors))
    left_bases = {
        name: _projector_basis_v233(
            left_symmetry_contract.projectors[name],
            tolerance=policy.projector_tolerance,
        )
        for name in names
    }
    right_bases = {
        name: _projector_basis_v233(
            right_symmetry_contract.projectors[name],
            tolerance=policy.projector_tolerance,
        )
        for name in names
    }
    blocks = []
    dimensions_match = True
    for name in names:
        left_basis = left_bases[name]
        right_basis = right_bases[name]
        dimensions_match = bool(
            dimensions_match and left_basis.shape[1] == right_basis.shape[1]
        )
        assigned = left_basis.conj().T @ overlap @ right_basis
        assigned_singular_values = np.linalg.svd(assigned, compute_uv=False)
        minimum = float(np.min(assigned_singular_values))
        competing = []
        for other in names:
            if other == name:
                continue
            block = left_basis.conj().T @ overlap @ right_bases[other]
            competing.append(float(np.linalg.norm(block, ord=2)))
        maximum_leakage = max(competing, default=0.0)
        blocks.append(
            ManifoldBlockTransportV233(
                name=name,
                dimension=int(left_basis.shape[1]),
                singular_values=np.asarray(assigned_singular_values, dtype=float),
                minimum_singular_value=minimum,
                maximum_competing_leakage=maximum_leakage,
                assignment_margin=float(minimum - maximum_leakage),
            )
        )

    J_left = np.asarray(
        left_symmetry_contract.time_reversal_matrix, dtype=complex
    )
    J_right = np.asarray(
        right_symmetry_contract.time_reversal_matrix, dtype=complex
    )
    time_reversal_residual = float(
        np.linalg.norm(
            overlap - J_left @ overlap.conj() @ J_right.conj().T,
            ord="fro",
        )
    )
    checks = {
        "left_complete_symmetry_contract": bool(left_audit.passed),
        "right_complete_symmetry_contract": bool(right_audit.passed),
        "same_complete_manifold_dimensions": bool(dimensions_match),
        "assigned_manifold_retention": bool(
            all(
                item.minimum_singular_value
                >= policy.minimum_assigned_singular_value
                for item in blocks
            )
        ),
        "competing_manifold_leakage": bool(
            all(
                item.maximum_competing_leakage
                <= policy.maximum_competing_leakage
                for item in blocks
            )
        ),
        "manifold_assignment_margin": bool(
            all(
                item.assignment_margin >= policy.minimum_assignment_margin
                for item in blocks
            )
        ),
        "time_reversal_covariance": bool(
            time_reversal_residual <= policy.time_reversal_tolerance
        ),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return CompleteManifoldTransportAuditV233(
        electron_parity=left_symmetry_contract.electron_parity,
        manifold_blocks=tuple(blocks),
        time_reversal_covariance_residual=time_reversal_residual,
        checks=checks,
        passed=bool(all(checks.values())),
        policy=policy.as_dict(),
    )


def require_complete_manifold_transport_v233(*args, **kwargs):
    report = audit_complete_manifold_transport_v233(*args, **kwargs)
    if not report.passed:
        failed = ", ".join(
            name for name, passed in report.checks.items() if not passed
        )
        raise ValueError("complete-manifold transport failed: " + failed)
    return report
