"""SOC-neutral electronic operator and provenance contract for v0.21.3.

This module introduces no physical spin-orbit Hamiltonian.  It freezes the interface
through which a later analytic or molecular backend may provide one.
"""

from dataclasses import dataclass, field, asdict
import hashlib
import json
import numpy as np

from .electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from .matrix_invariants_v213 import antihermiticity_residual_v213


HARTREE_PER_WAVENUMBER_V213 = 1.0 / 219474.63136320


def wavenumber_to_hartree_v213(value):
    return np.asarray(value) * HARTREE_PER_WAVENUMBER_V213


def hartree_to_wavenumber_v213(value):
    return np.asarray(value) / HARTREE_PER_WAVENUMBER_V213


def _canonical_json_value(value):
    if isinstance(value, np.generic):
        return _canonical_json_value(value.item())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("provenance cannot contain non-finite complex values.")
        return {"__complex__": [float(value.real), float(value.imag)]}
    if isinstance(value, np.ndarray):
        return _canonical_json_value(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("provenance dictionary keys must be strings.")
        return {
            key: _canonical_json_value(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_json_value(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("provenance cannot contain non-finite floating-point values.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported provenance value {type(value).__name__}.")


@dataclass(frozen=True)
class ElectronicStateDescriptorV213:
    label: str
    source_root: str | None = None
    multiplicity: int | None = None
    component: str | None = None
    charge: int | None = None

    def validate(self):
        if not str(self.label).strip():
            raise ValueError("every electronic state requires a nonempty label.")
        if self.multiplicity is not None:
            multiplicity = int(self.multiplicity)
            if multiplicity < 1 or multiplicity != self.multiplicity:
                raise ValueError("multiplicity must be a positive integer.")
        if self.source_root is not None and not str(self.source_root).strip():
            raise ValueError("source_root cannot be empty when supplied.")
        if self.component is not None and not str(self.component).strip():
            raise ValueError("component cannot be empty when supplied.")
        if self.charge is not None:
            charge = int(self.charge)
            if charge != self.charge:
                raise ValueError("charge must be an integer.")
        return self

    def as_dict(self):
        self.validate()
        return asdict(self)


@dataclass(frozen=True)
class ElectronicModelSpaceV213:
    name: str
    representation: str
    states: tuple[ElectronicStateDescriptorV213, ...]
    fixed_dimension: bool = True
    full_electronic_blocks: bool = True
    complete_multiplets: bool = False
    energy_unit: str = "hartree"
    coordinate_unit: str = "bohr"

    def __post_init__(self):
        object.__setattr__(self, "states", tuple(self.states))

    @property
    def nstate(self):
        return len(self.states)

    @property
    def fixed_frame(self):
        return self.representation in {"fixed_general", "fixed_spin_diabatic"}

    def validate(self):
        allowed = {
            "fixed_general",
            "fixed_spin_diabatic",
            "local_general",
            "total_hamiltonian_diagonal",
        }
        if self.representation not in allowed:
            raise ValueError(f"unsupported electronic representation {self.representation!r}.")
        if not str(self.name).strip():
            raise ValueError("model-space name cannot be empty.")
        for name, value in (
            ("fixed_dimension", self.fixed_dimension),
            ("full_electronic_blocks", self.full_electronic_blocks),
            ("complete_multiplets", self.complete_multiplets),
        ):
            if not isinstance(value, (bool, np.bool_)):
                raise ValueError(f"{name} must be Boolean.")
        if not self.fixed_dimension:
            raise ValueError("v0.21.3 requires a fixed electronic model-space dimension.")
        if not self.full_electronic_blocks:
            raise ValueError("v0.21.3 never prunes within an electronic block.")
        if self.energy_unit != "hartree" or self.coordinate_unit != "bohr":
            raise ValueError("internal operator data must use hartree and bohr.")
        if not self.states:
            raise ValueError("electronic model space cannot be empty.")
        for state in self.states:
            state.validate()
        labels = [state.label for state in self.states]
        if len(set(labels)) != len(labels):
            raise ValueError("electronic-state labels must be unique.")

        if self.complete_multiplets:
            groups = {}
            for state in self.states:
                if state.source_root is None or state.multiplicity is None:
                    raise ValueError(
                        "complete-multiplet validation requires source_root and multiplicity."
                    )
                key = (state.source_root, int(state.multiplicity), state.charge)
                groups.setdefault(key, []).append(state)
            for (root, multiplicity, _), group in groups.items():
                if len(group) != multiplicity:
                    raise ValueError(
                        f"root {root!r} has {len(group)} components; "
                        f"multiplicity {multiplicity} requires {multiplicity}."
                    )
                components = [state.component for state in group]
                if any(component is None for component in components):
                    raise ValueError("complete multiplets require explicit component labels.")
                if len(set(components)) != len(components):
                    raise ValueError("multiplet-component labels must be unique per root.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "name": self.name,
            "representation": self.representation,
            "states": [state.as_dict() for state in self.states],
            "fixed_dimension": self.fixed_dimension,
            "full_electronic_blocks": self.full_electronic_blocks,
            "complete_multiplets": self.complete_multiplets,
            "energy_unit": self.energy_unit,
            "coordinate_unit": self.coordinate_unit,
        }


@dataclass(frozen=True)
class ElectronicOperatorProvenanceV213:
    model_name: str
    model_version: str
    model_space: ElectronicModelSpaceV213
    spin_free_method: str
    soc_enabled: bool = False
    soc_method: str = "none"
    scalar_relativistic_method: str = "none"
    derivative_method: str = "analytic"
    parameters: dict = field(default_factory=dict)

    def validate(self):
        self.model_space.validate()
        if not isinstance(self.soc_enabled, (bool, np.bool_)):
            raise ValueError("soc_enabled must be Boolean.")
        if not isinstance(self.parameters, dict):
            raise ValueError("operator provenance parameters must be a dictionary.")
        for name, value in (
            ("model_name", self.model_name),
            ("model_version", self.model_version),
            ("spin_free_method", self.spin_free_method),
            ("soc_method", self.soc_method),
            ("scalar_relativistic_method", self.scalar_relativistic_method),
            ("derivative_method", self.derivative_method),
        ):
            if not str(value).strip():
                raise ValueError(f"{name} cannot be empty.")
        if self.soc_enabled and self.soc_method == "none":
            raise ValueError("soc_enabled=True requires an explicit soc_method.")
        if not self.soc_enabled and self.soc_method != "none":
            raise ValueError("soc_enabled=False requires soc_method='none'.")
        _canonical_json_value(self.parameters)
        return self

    def as_dict(self):
        self.validate()
        return {
            "model_name": self.model_name,
            "model_version": self.model_version,
            "model_space": self.model_space.as_dict(),
            "spin_free_method": self.spin_free_method,
            "soc_enabled": bool(self.soc_enabled),
            "soc_method": self.soc_method,
            "scalar_relativistic_method": self.scalar_relativistic_method,
            "derivative_method": self.derivative_method,
            "parameters": _canonical_json_value(self.parameters),
            "derivative_semantics": (
                "K[a]=<Phi|partial_a(H_spin_free+H_SOC)|Phi>; "
                "K is not the naive derivative of a moving-frame matrix"
            ),
        }

    def fingerprint(self):
        payload = json.dumps(
            _canonical_json_value(self.as_dict()),
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


def compose_electronic_operator_v213(
    *,
    q,
    H_spin_free,
    dH_spin_free_dq,
    connection_q,
    mass_matrix_q_au,
    provenance,
    H_soc=None,
    dH_soc_dq=None,
):
    """Compose total H and physical derivative operators without choosing a SOC model."""
    provenance = provenance.validate()
    H0 = np.asarray(H_spin_free, dtype=complex)
    K0 = np.asarray(dH_spin_free_dq, dtype=complex)
    Hso = np.zeros_like(H0) if H_soc is None else np.asarray(H_soc, dtype=complex)
    Kso = np.zeros_like(K0) if dH_soc_dq is None else np.asarray(dH_soc_dq, dtype=complex)
    if Hso.shape != H0.shape or Kso.shape != K0.shape:
        raise ValueError("spin-free and SOC operator terms must have matching shapes.")
    soc_signal = max(np.linalg.norm(Hso, ord="fro"), np.linalg.norm(Kso.ravel()))
    if soc_signal > 0.0 and not provenance.soc_enabled:
        raise ValueError("nonzero SOC terms require soc_enabled provenance.")
    if provenance.model_space.nstate != H0.shape[0]:
        raise ValueError("operator dimension does not match the declared model space.")

    metadata = {
        "v213_electronic_contract": provenance.as_dict(),
        "v213_provenance_fingerprint": provenance.fingerprint(),
        "operator_decomposition": {
            "H": "H_spin_free + H_SOC",
            "K": "dH_spin_free_dq + dH_SOC_dq",
            "soc_enabled": bool(provenance.soc_enabled),
        },
    }
    point = ElectronicOperatorPointV21(
        q=np.asarray(q),
        H=H0 + Hso,
        dH_dq=K0 + Kso,
        connection_q=np.asarray(connection_q, dtype=complex),
        mass_matrix_q_au=np.asarray(mass_matrix_q_au),
        metadata=metadata,
    ).validate()
    validate_electronic_contract_v213(point, provenance)
    return point


def validate_electronic_contract_v213(point, provenance, *, tolerance=1.0e-12):
    point = point.validate(atol=tolerance)
    provenance = provenance.validate()
    if point.nstate != provenance.model_space.nstate:
        raise ValueError("electronic dimension differs from the declared model space.")
    if provenance.model_space.fixed_frame:
        residuals = [
            antihermiticity_residual_v213(point.connection_q[a])
            for a in range(point.nq)
        ]
        magnitude = max(
            (float(np.linalg.norm(point.connection_q[a], ord="fro")) for a in range(point.nq)),
            default=0.0,
        )
        if magnitude > tolerance or max(residuals, default=0.0) > tolerance:
            raise ValueError("a fixed electronic frame requires connection_q=0.")
    metadata = dict(point.metadata)
    fingerprint = metadata.get("v213_provenance_fingerprint")
    if fingerprint is not None and fingerprint != provenance.fingerprint():
        raise ValueError("operator metadata and requested provenance disagree.")
    return point


class ContractedElectronicOperatorProviderV213:
    """Validate dimension, units, derivative semantics, and provenance on every call."""

    def __init__(self, base_provider, provenance):
        self.base_provider = base_provider
        self.provenance = provenance.validate()
        self._fingerprint = self.provenance.fingerprint()

    def evaluate_snapshot(self, q):
        base = self.base_provider.evaluate_snapshot(q)
        metadata = {
            **dict(base.point.metadata),
            "v213_electronic_contract": self.provenance.as_dict(),
            "v213_provenance_fingerprint": self._fingerprint,
        }
        point = ElectronicOperatorPointV21(
            q=np.asarray(base.point.q, dtype=float).copy(),
            H=np.asarray(base.point.H, dtype=complex).copy(),
            dH_dq=np.asarray(base.point.dH_dq, dtype=complex).copy(),
            connection_q=np.asarray(base.point.connection_q, dtype=complex).copy(),
            mass_matrix_q_au=np.asarray(base.point.mass_matrix_q_au, dtype=float).copy(),
            metadata=metadata,
        )
        validate_electronic_contract_v213(point, self.provenance)
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=(
                None
                if base.state_vectors is None
                else np.asarray(base.state_vectors, dtype=complex).copy()
            ),
            wavefunction_snapshot=base.wavefunction_snapshot,
            parent_snapshot=base,
            frame_from_parent=np.eye(point.nstate, dtype=complex),
            metadata={
                **dict(base.metadata),
                "provider": "ContractedElectronicOperatorProviderV213",
                "provenance_fingerprint": self._fingerprint,
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        if left.state_vectors is not None and right.state_vectors is not None:
            return left.state_vectors.conj().T @ right.state_vectors
        return self.base_provider.snapshot_overlap(
            left.parent_snapshot,
            right.parent_snapshot,
        )

    def diagnostics_dict(self):
        base = (
            self.base_provider.diagnostics_dict()
            if hasattr(self.base_provider, "diagnostics_dict")
            else {}
        )
        return {
            "v213_contract": self.provenance.as_dict(),
            "provenance_fingerprint": self._fingerprint,
            "base": base,
        }
