"""Structural and runtime hardening for method-specific PySCF SOC engines.

No physical SOC implementation is supplied here.  A concrete engine remains
responsible for every electronic-structure quantity and raw-artifact parser.
"""

from dataclasses import asdict, dataclass

import numpy as np

from .molecular_soc_dossier_v231 import BackendRuntimeAttestationV231
from .molecular_soc_runtime_v232 import (
    BackendArtifactValidationProofV232,
    BackendMethodIdentityV232,
    convergence_from_snapshot_v232,
)
from .pyscf_nac_convention_v232 import PYSCF_NAC_EMPIRICAL_MAPPING_V232
from .pyscf_runtime_v232 import require_pyscf_runtime_v232


PYSCF_NAC_CONVENTION_V232 = PYSCF_NAC_EMPIRICAL_MAPPING_V232
_REQUIRED_CALLABLES_V232 = (
    "components",
    "evaluate_snapshot",
    "snapshot_overlap",
    "write_raw_artifacts",
    "validate_raw_artifacts_v232",
)


def _native_bool_v232(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


@dataclass(frozen=True)
class PySCFMethodSpecificCapabilitiesV232:
    state_interaction_soc: bool
    physical_soc_derivatives: bool
    analytic_spin_free_gradients: bool
    derivative_connections: bool
    many_electron_overlaps: bool
    raw_artifact_parser: bool
    fresh_execution: bool

    def validate(self):
        for name in asdict(self):
            _native_bool_v232(name, getattr(self, name))
        return self

    @property
    def complete(self):
        self.validate()
        return bool(all(asdict(self).values()))


def validate_pyscf_engine_contract_v232(
    method_specific_engine,
    *,
    runtime_version,
    expected_identity,
):
    """Validate structure and exact identity before consulting capability flags."""
    missing = []
    noncallable = []
    for name in _REQUIRED_CALLABLES_V232:
        if not hasattr(method_specific_engine, name):
            missing.append(name)
        elif not callable(getattr(method_specific_engine, name)):
            noncallable.append(name)
    if missing or noncallable:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if noncallable:
            details.append("noncallable " + ", ".join(noncallable))
        raise TypeError("method-specific PySCF SOC engine is incomplete: " + "; ".join(details))

    if type(expected_identity) is not BackendMethodIdentityV232:
        raise TypeError("expected_identity must be BackendMethodIdentityV232.")
    expected_identity.validate()
    engine_identity = getattr(method_specific_engine, "method_identity", None)
    if type(engine_identity) is not BackendMethodIdentityV232:
        raise TypeError("engine method_identity must be BackendMethodIdentityV232.")
    engine_identity.validate()
    if engine_identity != expected_identity:
        raise ValueError("engine method identity differs from the trusted identity.")
    if engine_identity.backend_name != "PySCF":
        raise ValueError("PySCF adapter requires exact backend_name='PySCF'.")
    if engine_identity.source_kind != "live_ab_initio":
        raise ValueError("PySCF live adapter requires source_kind='live_ab_initio'.")
    if engine_identity.backend_version != runtime_version:
        raise ValueError("engine version differs from the imported PySCF runtime.")
    if engine_identity.nac_convention != PYSCF_NAC_CONVENTION_V232:
        raise ValueError("engine uses the wrong PySCF NAC orientation.")

    capabilities = getattr(method_specific_engine, "capabilities", None)
    if type(capabilities) is not PySCFMethodSpecificCapabilitiesV232:
        raise TypeError(
            "engine capabilities must be PySCFMethodSpecificCapabilitiesV232."
        )
    capabilities.validate()
    if not capabilities.complete:
        raise RuntimeError("method-specific PySCF SOC capabilities are incomplete.")
    return capabilities


def _installed_pyscf_version_v232():
    probe = require_pyscf_runtime_v232()
    return str(probe.module_version)


class PySCFMethodSpecificSOCAdapterV232:
    """Expose a structurally complete engine with exact, trusted method identity."""

    def __init__(self, method_specific_engine, *, expected_identity):
        runtime_version = _installed_pyscf_version_v232()
        capabilities = validate_pyscf_engine_contract_v232(
            method_specific_engine,
            runtime_version=runtime_version,
            expected_identity=expected_identity,
        )
        self.engine = method_specific_engine
        self.method_identity = expected_identity
        self.capabilities = capabilities
        self.runtime_version = runtime_version

    @property
    def adapter_name(self):
        return self.method_identity.adapter_name

    @property
    def adapter_version(self):
        return self.method_identity.adapter_version

    def components(self, q):
        return self.engine.components(q).validate()

    def evaluate_snapshot(self, q):
        snapshot = self.engine.evaluate_snapshot(q).validate()
        convergence = convergence_from_snapshot_v232(snapshot)
        if not convergence.complete:
            failed = ", ".join(
                name
                for name, value in asdict(convergence).items()
                if name != "vocabulary" and not value
            )
            raise RuntimeError(
                "live PySCF snapshot has unconverged required stages: " + failed
            )
        return snapshot

    def snapshot_overlap(self, left, right):
        return self.engine.snapshot_overlap(left, right)

    def write_raw_artifacts(self, *args, **kwargs):
        return self.engine.write_raw_artifacts(*args, **kwargs)

    def validate_raw_artifacts_v232(self, *args, **kwargs):
        proof = self.engine.validate_raw_artifacts_v232(*args, **kwargs)
        if type(proof) is not BackendArtifactValidationProofV232:
            raise TypeError(
                "PySCF raw-artifact parser must return "
                "BackendArtifactValidationProofV232, not a Boolean."
            )
        return proof.validate()

    def runtime_attestation(self, *, environment_sha256, runtime_probe_artifact):
        """Build the legacy dossier attestation; v0.23.2 proof adds exact binding."""
        return BackendRuntimeAttestationV231(
            runtime_name=self.method_identity.backend_name,
            runtime_version=self.runtime_version,
            adapter_name=self.method_identity.adapter_name,
            adapter_version=self.method_identity.adapter_version,
            environment_sha256=environment_sha256,
            runtime_probe_artifact=runtime_probe_artifact,
            runtime_imported=True,
            method_specific_soc_implemented=self.capabilities.state_interaction_soc,
            soc_derivatives_implemented=self.capabilities.physical_soc_derivatives,
            wavefunction_overlaps_implemented=self.capabilities.many_electron_overlaps,
            artifact_parser_validated=self.capabilities.raw_artifact_parser,
            fresh_execution_observed=self.capabilities.fresh_execution,
        ).validate()
