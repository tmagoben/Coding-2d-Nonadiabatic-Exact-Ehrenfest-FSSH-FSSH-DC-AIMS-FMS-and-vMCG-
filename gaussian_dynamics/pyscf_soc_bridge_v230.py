"""Fail-closed boundary for a future live PySCF molecular SOC implementation."""

from dataclasses import asdict, dataclass
import importlib.util

from .molecular_soc_contract_v230 import (
    molecular_soc_contract_from_provider_v230,
)


@dataclass(frozen=True)
class PySCFSOCRuntimeProbeV230:
    installed: bool
    version: str | None
    live_soc_adapter_validated: bool
    note: str

    def as_dict(self):
        return asdict(self)


def probe_pyscf_soc_runtime_v230():
    """Report availability without inferring unsupported SOC capabilities."""
    if importlib.util.find_spec("pyscf") is None:
        return PySCFSOCRuntimeProbeV230(
            installed=False,
            version=None,
            live_soc_adapter_validated=False,
            note="PySCF is not installed; no live molecular SOC calculation was run.",
        )
    import pyscf

    return PySCFSOCRuntimeProbeV230(
        installed=True,
        version=str(getattr(pyscf, "__version__", "unknown")),
        live_soc_adapter_validated=False,
        note=(
            "PySCF import succeeded, but importability alone does not validate the "
            "state-resolved SOC derivative and overlap contract."
        ),
    )


def require_pyscf_soc_runtime_v230():
    probe = probe_pyscf_soc_runtime_v230()
    if not probe.installed:
        raise ImportError(
            "PySCF is required for a live molecular SOC bridge; no live backend "
            "is admitted by the deterministic replay fixtures."
        )
    return probe


class PySCFMolecularSOCBridgeV230:
    """Validate and expose an injected live PySCF SOC implementation.

    The framework deliberately does not guess how a particular PySCF method constructs
    state-interaction SOC matrices or their physical nuclear derivatives.  A method-
    specific implementation must expose the complete framework provider protocol and
    is then checked here before use.
    """

    def __init__(self, method_specific_provider):
        probe = require_pyscf_soc_runtime_v230()
        contract = molecular_soc_contract_from_provider_v230(method_specific_provider)
        identity = contract.identity
        if identity.backend_name.strip().lower() != "pyscf":
            raise ValueError("PySCF bridge requires backend_name='PySCF'.")
        if identity.source_kind != "live_ab_initio":
            raise ValueError("PySCF bridge requires source_kind='live_ab_initio'.")
        if identity.backend_version != probe.version:
            raise ValueError("declared PySCF version differs from the imported runtime.")
        required = (
            "components",
            "evaluate_snapshot",
            "snapshot_overlap",
            "provenance",
            "soc_symmetry_contract",
            "molecular_soc_contract",
        )
        missing = [name for name in required if not hasattr(method_specific_provider, name)]
        if missing:
            raise TypeError(
                "method-specific PySCF SOC provider is incomplete: " + ", ".join(missing)
            )
        self.base_provider = method_specific_provider
        self.provenance = method_specific_provider.provenance
        self._contract = contract

    @property
    def molecular_soc_contract(self):
        return self._contract

    @property
    def soc_symmetry_contract(self):
        return self.base_provider.soc_symmetry_contract

    @property
    def time_reversal_matrix(self):
        return self.base_provider.time_reversal_matrix

    @property
    def projectors(self):
        return self.base_provider.projectors

    def components(self, q):
        return self.base_provider.components(q).validate()

    def evaluate_snapshot(self, q):
        snapshot = self.base_provider.evaluate_snapshot(q).validate()
        metadata = {**dict(snapshot.point.metadata), **dict(snapshot.metadata)}
        required_convergence = (
            "scf_converged",
            "correlated_converged",
            "soc_converged",
        )
        if not all(metadata.get(name) is True for name in required_convergence):
            raise RuntimeError(
                "live PySCF SOC snapshot lacks affirmative SCF, correlated, and SOC "
                "convergence evidence."
            )
        return snapshot

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        return self.base_provider.snapshot_overlap(left, right)
