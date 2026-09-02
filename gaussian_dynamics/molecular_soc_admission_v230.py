"""Admission audit separating protocol validity from real molecular evidence."""

from dataclasses import asdict, dataclass
import numpy as np

from .molecular_soc_contract_v230 import (
    molecular_soc_contract_from_provider_v230,
    require_trajectory_ready_molecular_soc_v230,
)
from .physical_soc_validation_v220 import (
    PhysicalSOCAuditSettingsV220,
    audit_physical_soc_provider_v220,
)
from .soc_admission_v221 import (
    audit_soc_symmetry_contract_v221,
    soc_symmetry_contract_from_provider_v221,
)


@dataclass(frozen=True)
class MolecularSOCAdmissionSettingsV230:
    physical_soc_settings: PhysicalSOCAuditSettingsV220 = PhysicalSOCAuditSettingsV220()
    symmetry_tolerance: float = 1.0e-12

    def validate(self):
        self.physical_soc_settings.validate()
        if not np.isfinite(self.symmetry_tolerance) or self.symmetry_tolerance <= 0.0:
            raise ValueError("molecular SOC symmetry tolerance must be finite and positive.")
        return self


@dataclass(frozen=True)
class MolecularSOCAdmissionAuditV230:
    q: np.ndarray
    backend_name: str
    backend_version: str
    source_kind: str
    capability_tier: str
    protocol_checks: dict
    real_admission_checks: dict
    protocol_passed: bool
    real_backend_admitted: bool
    passed: bool
    physical_soc_report: object
    symmetry_report: object
    provenance_fingerprint: str
    molecular_contract_fingerprint: str
    thresholds: dict

    def as_dict(self):
        return {
            "q": np.asarray(self.q, dtype=float).tolist(),
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "source_kind": self.source_kind,
            "capability_tier": self.capability_tier,
            "protocol_checks": dict(self.protocol_checks),
            "real_admission_checks": dict(self.real_admission_checks),
            "protocol_passed": bool(self.protocol_passed),
            "real_backend_admitted": bool(self.real_backend_admitted),
            "passed": bool(self.passed),
            "physical_soc_report": self.physical_soc_report.as_dict(),
            "symmetry_report": self.symmetry_report.as_dict(),
            "provenance_fingerprint": self.provenance_fingerprint,
            "molecular_contract_fingerprint": self.molecular_contract_fingerprint,
            "thresholds": dict(self.thresholds),
        }


def audit_molecular_soc_provider_v230(
    provider,
    q,
    *,
    require_real_backend=True,
    settings=MolecularSOCAdmissionSettingsV230(),
):
    """Audit a replay or live provider without conflating fixtures with ab initio data."""
    settings = settings.validate()
    if not isinstance(require_real_backend, (bool, np.bool_)):
        raise ValueError("require_real_backend must be Boolean.")
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or len(q) < 1 or not np.all(np.isfinite(q)):
        raise ValueError("molecular SOC admission requires a finite coordinate vector.")
    contract = molecular_soc_contract_from_provider_v230(provider)
    trajectory_contract = require_trajectory_ready_molecular_soc_v230(provider)
    symmetry_contract = soc_symmetry_contract_from_provider_v221(provider)
    provenance = provider.provenance.validate()
    symmetry = audit_soc_symmetry_contract_v221(
        provenance.model_space,
        symmetry_contract,
        provenance=provenance,
        fermionic=(contract.identity.electron_parity == "odd"),
        tolerance=settings.symmetry_tolerance,
    )
    physical = audit_physical_soc_provider_v220(
        provider,
        q,
        fermionic=(contract.identity.electron_parity == "odd"),
        settings=settings.physical_soc_settings,
    )
    snapshot = provider.evaluate_snapshot(q).validate()
    parameters = provenance.parameters
    stored_contract = parameters.get("v230_molecular_soc_contract")
    stored_fingerprint = parameters.get("v230_molecular_soc_contract_fingerprint")
    replay_integrity = bool(
        not hasattr(provider, "replay_fingerprint")
        or (
            isinstance(provider.replay_fingerprint, str)
            and len(provider.replay_fingerprint) == 64
            and snapshot.point.metadata.get("v230_replay_dataset_fingerprint")
            == provider.replay_fingerprint
        )
    )
    emitted_convergence = snapshot.point.metadata.get(
        "v230_electronic_converged",
        contract.all_electronic_calculations_converged,
    )
    protocol_checks = {
        "static_SOC_capability": bool(contract.capabilities.static_soc),
        "trajectory_capabilities": bool(trajectory_contract.capabilities.trajectory_ready),
        "all_electronic_calculations_converged": bool(
            contract.all_electronic_calculations_converged
            and bool(emitted_convergence)
        ),
        "single_electron_parity_and_charge": bool(symmetry.passed),
        "physical_SOC_contract": bool(physical.passed),
        "component_resolved_derivatives": bool(
            physical.checks["spin_free_component_derivatives"]
            and physical.checks["SOC_component_derivatives"]
        ),
        "cross_geometry_differentials": bool(
            physical.checks["cross_geometry_differentials"]
        ),
        "molecular_contract_provenance": bool(
            stored_contract == contract.as_dict()
            and stored_fingerprint == contract.fingerprint()
        ),
        "replay_integrity": replay_integrity,
    }
    evidence = contract.evidence
    real_checks = {
        **protocol_checks,
        "real_ab_initio_source": bool(contract.identity.real_ab_initio_source),
        "independent_reference_evidence": bool(
            evidence.independent_reference_validated
        ),
        "basis_convergence_evidence": bool(evidence.basis_converged),
        "method_convergence_evidence": bool(evidence.method_converged),
        "translation_rotation_invariance": bool(
            evidence.frame_invariance_validated
        ),
        "state_tracking_quality": bool(evidence.state_tracking_validated),
        "traceable_nuclear_identity": bool(
            contract.identity.traceable_nuclear_identity
        ),
        "contract_real_admission_ready": bool(contract.real_backend_admission_ready),
    }
    protocol_passed = bool(all(protocol_checks.values()))
    real_admitted = bool(all(real_checks.values()))
    return MolecularSOCAdmissionAuditV230(
        q=q.copy(),
        backend_name=contract.identity.backend_name,
        backend_version=contract.identity.backend_version,
        source_kind=contract.identity.source_kind,
        capability_tier=contract.capabilities.tier,
        protocol_checks={name: bool(value) for name, value in protocol_checks.items()},
        real_admission_checks={name: bool(value) for name, value in real_checks.items()},
        protocol_passed=protocol_passed,
        real_backend_admitted=real_admitted,
        passed=real_admitted if require_real_backend else protocol_passed,
        physical_soc_report=physical,
        symmetry_report=symmetry,
        provenance_fingerprint=provenance.fingerprint(),
        molecular_contract_fingerprint=contract.fingerprint(),
        thresholds=asdict(settings),
    )


def require_molecular_soc_protocol_v230(*args, **kwargs):
    kwargs["require_real_backend"] = False
    report = audit_molecular_soc_provider_v230(*args, **kwargs)
    if not report.protocol_passed:
        failed = ", ".join(
            name for name, value in report.protocol_checks.items() if not value
        )
        raise ValueError(f"molecular SOC protocol failed: {failed}.")
    return report


def require_real_molecular_soc_backend_v230(*args, **kwargs):
    kwargs["require_real_backend"] = True
    report = audit_molecular_soc_provider_v230(*args, **kwargs)
    if not report.real_backend_admitted:
        failed = ", ".join(
            name for name, value in report.real_admission_checks.items() if not value
        )
        raise ValueError(f"real molecular SOC backend admission failed: {failed}.")
    return report
