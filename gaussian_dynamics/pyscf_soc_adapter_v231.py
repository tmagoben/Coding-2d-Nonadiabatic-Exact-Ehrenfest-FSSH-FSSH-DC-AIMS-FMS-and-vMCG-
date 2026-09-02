"""Fail-closed method-specific PySCF SOC adapter boundary for v0.23.1.

PySCF importability and its documented NAC API do not by themselves establish a
state-interaction SOC method, physical SOC derivatives, or many-electron overlap
tracking.  Those capabilities must be supplied and validated by a concrete engine.
"""

from dataclasses import asdict, dataclass
import importlib.util
import numpy as np

from .molecular_soc_dossier_v231 import BackendRuntimeAttestationV231


PYSCF_NAC_CONVENTION_V231 = "state=(ket,bra) returns <bra|d ket/dR>"


def _native_bool_v231(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


@dataclass(frozen=True)
class PySCFMethodSpecificCapabilitiesV231:
    state_interaction_soc: bool
    physical_soc_derivatives: bool
    analytic_spin_free_gradients: bool
    derivative_connections: bool
    many_electron_overlaps: bool
    raw_artifact_parser: bool
    fresh_execution: bool

    def validate(self):
        for name in (
            "state_interaction_soc",
            "physical_soc_derivatives",
            "analytic_spin_free_gradients",
            "derivative_connections",
            "many_electron_overlaps",
            "raw_artifact_parser",
            "fresh_execution",
        ):
            _native_bool_v231(name, getattr(self, name))
        return self

    @property
    def complete(self):
        self.validate()
        return bool(all(asdict(self).values()))


@dataclass(frozen=True)
class PySCFSOCAdapterProbeV231:
    installed: bool
    version: str | None
    core_nac_api_available: bool
    method_specific_adapter_supplied: bool
    live_admission_ready: bool
    nac_convention: str
    note: str

    def as_dict(self):
        return asdict(self)


def probe_pyscf_soc_adapter_v231(method_specific_engine=None):
    if importlib.util.find_spec("pyscf") is None:
        return PySCFSOCAdapterProbeV231(
            installed=False,
            version=None,
            core_nac_api_available=False,
            method_specific_adapter_supplied=method_specific_engine is not None,
            live_admission_ready=False,
            nac_convention=PYSCF_NAC_CONVENTION_V231,
            note=(
                "PySCF is not installed; no live SOC calculation or method-specific "
                "artifact validation was run."
            ),
        )
    import pyscf

    try:
        from pyscf.nac.sacasscf import NonAdiabaticCouplings  # noqa: F401

        nac_available = True
    except (ImportError, AttributeError):
        nac_available = False
    complete = bool(
        method_specific_engine is not None
        and isinstance(
            getattr(method_specific_engine, "capabilities", None),
            PySCFMethodSpecificCapabilitiesV231,
        )
        and method_specific_engine.capabilities.complete
        and getattr(method_specific_engine, "nac_convention", None)
        == PYSCF_NAC_CONVENTION_V231
    )
    return PySCFSOCAdapterProbeV231(
        installed=True,
        version=str(getattr(pyscf, "__version__", "unknown")),
        core_nac_api_available=nac_available,
        method_specific_adapter_supplied=method_specific_engine is not None,
        live_admission_ready=bool(nac_available and complete),
        nac_convention=PYSCF_NAC_CONVENTION_V231,
        note=(
            "PySCF and its SA-CASSCF NAC API are available; live admission additionally "
            "requires a complete method-specific SOC engine."
            if nac_available
            else "PySCF is importable but the required SA-CASSCF NAC API is unavailable."
        ),
    )


def require_pyscf_soc_adapter_v231(method_specific_engine=None):
    probe = probe_pyscf_soc_adapter_v231(method_specific_engine)
    if not probe.installed:
        raise ImportError("PySCF is not installed; live v0.23.1 SOC admission is unavailable.")
    if not probe.core_nac_api_available:
        raise RuntimeError("the installed PySCF lacks the required SA-CASSCF NAC API.")
    if not probe.method_specific_adapter_supplied:
        raise RuntimeError(
            "a method-specific PySCF state-interaction SOC engine is required."
        )
    if not probe.live_admission_ready:
        raise RuntimeError("the method-specific PySCF SOC engine is incomplete.")
    return probe


class PySCFMethodSpecificSOCAdapterV231:
    """Expose only a fully declared, externally implemented PySCF SOC engine."""

    def __init__(self, method_specific_engine):
        probe = require_pyscf_soc_adapter_v231(method_specific_engine)
        required_attributes = (
            "backend_version",
            "adapter_name",
            "adapter_version",
            "nac_convention",
            "components",
            "evaluate_snapshot",
            "snapshot_overlap",
            "write_raw_artifacts",
            "validate_raw_artifacts",
        )
        missing = [
            name for name in required_attributes if not hasattr(method_specific_engine, name)
        ]
        if missing:
            raise TypeError(
                "method-specific PySCF SOC engine is incomplete: " + ", ".join(missing)
            )
        if str(method_specific_engine.backend_version) != probe.version:
            raise ValueError("method-specific engine version differs from imported PySCF.")
        if method_specific_engine.nac_convention != PYSCF_NAC_CONVENTION_V231:
            raise ValueError("method-specific engine uses the wrong PySCF NAC convention.")
        self.engine = method_specific_engine
        self.probe = probe
        self.capabilities = method_specific_engine.capabilities.validate()

    @property
    def adapter_name(self):
        return str(self.engine.adapter_name)

    @property
    def adapter_version(self):
        return str(self.engine.adapter_version)

    def components(self, q):
        return self.engine.components(q).validate()

    def evaluate_snapshot(self, q):
        snapshot = self.engine.evaluate_snapshot(q).validate()
        metadata = {**dict(snapshot.point.metadata), **dict(snapshot.metadata)}
        required = (
            "scf_converged",
            "correlated_converged",
            "soc_converged",
            "derivatives_converged",
            "overlaps_converged",
        )
        if not all(metadata.get(name) is True for name in required):
            raise RuntimeError(
                "live PySCF snapshot lacks affirmative convergence for every required stage."
            )
        return snapshot

    def snapshot_overlap(self, left, right):
        return self.engine.snapshot_overlap(left, right)

    def write_raw_artifacts(self, *args, **kwargs):
        return self.engine.write_raw_artifacts(*args, **kwargs)

    def validate_raw_artifacts(self, *args, **kwargs):
        result = self.engine.validate_raw_artifacts(*args, **kwargs)
        if result is not True:
            raise RuntimeError("method-specific PySCF raw-artifact validation failed.")
        return True

    def runtime_attestation(self, *, environment_sha256, runtime_probe_artifact):
        return BackendRuntimeAttestationV231(
            runtime_name="PySCF",
            runtime_version=self.probe.version,
            adapter_name=str(self.engine.adapter_name),
            adapter_version=str(self.engine.adapter_version),
            environment_sha256=environment_sha256,
            runtime_probe_artifact=runtime_probe_artifact,
            runtime_imported=True,
            method_specific_soc_implemented=self.capabilities.state_interaction_soc,
            soc_derivatives_implemented=self.capabilities.physical_soc_derivatives,
            wavefunction_overlaps_implemented=self.capabilities.many_electron_overlaps,
            artifact_parser_validated=self.capabilities.raw_artifact_parser,
            fresh_execution_observed=self.capabilities.fresh_execution,
        ).validate()
