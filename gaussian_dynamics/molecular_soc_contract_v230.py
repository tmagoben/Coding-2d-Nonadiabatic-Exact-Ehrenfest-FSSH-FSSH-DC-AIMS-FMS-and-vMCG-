"""Molecular spin-orbit backend admission contracts for v0.23.0.

The contract distinguishes a provider's declared numerical capabilities from evidence
that a real molecular backend has actually been validated.  Analytic fixtures may test
the protocol, but they can never satisfy the real-backend admission gate.
"""

from dataclasses import dataclass, field
import hashlib
import json
import numpy as np

from .electronic_contract_v213 import ElectronicOperatorProvenanceV213
from .soc_admission_v221 import SOCSymmetryContractV221


_SOURCE_KINDS_V230 = {
    "validation_fixture",
    "external_ab_initio_snapshot",
    "live_ab_initio",
}


def _native_bool_v230(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


def _required_text_v230(name, value):
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} cannot be empty.")
    return text


def _optional_nonnegative_v230(name, value):
    if value is None:
        return None
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative.")
    return result


def _canonical_v230(value):
    if isinstance(value, np.generic):
        return _canonical_v230(value.item())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("molecular SOC identity cannot contain non-finite values.")
        return {"__complex__": [float(value.real), float(value.imag)]}
    if isinstance(value, np.ndarray):
        return _canonical_v230(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("molecular SOC identity keys must be strings.")
        return {
            key: _canonical_v230(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v230(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("molecular SOC identity cannot contain non-finite values.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported molecular SOC identity value {type(value).__name__}.")


@dataclass(frozen=True)
class MolecularSOCCapabilitiesV230:
    static_soc: bool
    spin_free_derivatives: bool = False
    soc_derivatives: bool = False
    derivative_connections: bool = False
    cross_geometry_overlaps: bool = False
    deterministic_replay: bool = False
    analytic_soc_derivatives: bool = False

    def validate(self):
        values = {
            name: _native_bool_v230(name, getattr(self, name))
            for name in (
                "static_soc",
                "spin_free_derivatives",
                "soc_derivatives",
                "derivative_connections",
                "cross_geometry_overlaps",
                "deterministic_replay",
                "analytic_soc_derivatives",
            )
        }
        if not values["static_soc"]:
            raise ValueError("a molecular SOC contract must provide static SOC.")
        if values["analytic_soc_derivatives"] and not values["soc_derivatives"]:
            raise ValueError(
                "analytic_soc_derivatives requires the SOC-derivative capability."
            )
        return self

    @property
    def trajectory_ready(self):
        self.validate()
        return bool(
            self.static_soc
            and self.spin_free_derivatives
            and self.soc_derivatives
            and self.derivative_connections
            and self.cross_geometry_overlaps
        )

    @property
    def tier(self):
        return "trajectory_ready" if self.trajectory_ready else "static_soc"

    def as_dict(self):
        self.validate()
        return {
            "static_soc": bool(self.static_soc),
            "spin_free_derivatives": bool(self.spin_free_derivatives),
            "soc_derivatives": bool(self.soc_derivatives),
            "derivative_connections": bool(self.derivative_connections),
            "cross_geometry_overlaps": bool(self.cross_geometry_overlaps),
            "deterministic_replay": bool(self.deterministic_replay),
            "analytic_soc_derivatives": bool(self.analytic_soc_derivatives),
            "tier": self.tier,
        }


@dataclass(frozen=True)
class MolecularSOCBackendIdentityV230:
    backend_name: str
    backend_version: str
    source_kind: str
    electronic_method: str
    basis: str
    charge: int
    electron_count: int
    soc_operator: str
    scalar_relativistic_method: str
    derivative_method: str
    active_space: str = "none"
    molecule_name: str = "validation fixture"
    atom_symbols: tuple[str, ...] = ()
    isotope_masses_amu: tuple[float, ...] = ()
    reference_geometry_bohr: tuple[tuple[float, float, float], ...] = ()
    calculation_input_sha256: str | None = None
    environment_sha256: str | None = None
    geometry_unit: str = "bohr"
    energy_unit: str = "hartree"
    derivative_unit: str = "hartree/bohr"
    extra: dict = field(default_factory=dict)

    def validate(self):
        for name in (
            "backend_name",
            "backend_version",
            "electronic_method",
            "basis",
            "soc_operator",
            "scalar_relativistic_method",
            "derivative_method",
            "active_space",
            "molecule_name",
        ):
            _required_text_v230(name, getattr(self, name))
        if self.source_kind not in _SOURCE_KINDS_V230:
            raise ValueError(f"unsupported molecular SOC source_kind {self.source_kind!r}.")
        if self.geometry_unit != "bohr":
            raise ValueError("molecular SOC geometry data must use bohr.")
        if self.energy_unit != "hartree":
            raise ValueError("molecular SOC energies must use hartree.")
        if self.derivative_unit != "hartree/bohr":
            raise ValueError("molecular SOC derivatives must use hartree/bohr.")
        if int(self.charge) != self.charge:
            raise ValueError("molecular charge must be an integer.")
        if int(self.electron_count) != self.electron_count or int(self.electron_count) < 1:
            raise ValueError("electron_count must be a positive integer.")
        if not isinstance(self.extra, dict):
            raise ValueError("backend identity extra data must be a dictionary.")
        symbols = tuple(str(value).strip() for value in self.atom_symbols)
        masses = np.asarray(self.isotope_masses_amu, dtype=float)
        geometry = np.asarray(self.reference_geometry_bohr, dtype=float)
        if any(not value for value in symbols):
            raise ValueError("atom_symbols cannot contain empty labels.")
        if symbols:
            if masses.shape != (len(symbols),):
                raise ValueError("isotope_masses_amu must contain one mass per atom.")
            if geometry.shape != (len(symbols), 3):
                raise ValueError(
                    "reference_geometry_bohr must have shape (natom,3)."
                )
            if not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
                raise ValueError("all isotope masses must be finite and positive.")
            if not np.all(np.isfinite(geometry)):
                raise ValueError("reference molecular geometry must be finite.")
        elif masses.size or geometry.size:
            raise ValueError("nuclear masses or geometry require atom_symbols.")
        for name, fingerprint in (
            ("calculation_input_sha256", self.calculation_input_sha256),
            ("environment_sha256", self.environment_sha256),
        ):
            if fingerprint is not None:
                fingerprint = str(fingerprint)
                if len(fingerprint) != 64 or any(
                    character not in "0123456789abcdef" for character in fingerprint
                ):
                    raise ValueError(f"{name} must be a lowercase SHA-256 digest.")
        if self.real_ab_initio_source and not self.traceable_nuclear_identity:
            raise ValueError(
                "real molecular SOC sources require nuclear identity, reference "
                "geometry, and calculation/environment SHA-256 fingerprints."
            )
        _canonical_v230(self.extra)
        return self

    @property
    def real_ab_initio_source(self):
        return self.source_kind in {
            "external_ab_initio_snapshot",
            "live_ab_initio",
        }

    @property
    def electron_parity(self):
        return "odd" if int(self.electron_count) % 2 else "even"

    @property
    def traceable_nuclear_identity(self):
        return bool(
            self.atom_symbols
            and len(self.atom_symbols) == len(self.isotope_masses_amu)
            and len(self.atom_symbols) == len(self.reference_geometry_bohr)
            and self.calculation_input_sha256 is not None
            and self.environment_sha256 is not None
        )

    def as_dict(self):
        self.validate()
        return {
            "backend_name": self.backend_name,
            "backend_version": self.backend_version,
            "source_kind": self.source_kind,
            "electronic_method": self.electronic_method,
            "basis": self.basis,
            "charge": int(self.charge),
            "electron_count": int(self.electron_count),
            "electron_parity": self.electron_parity,
            "soc_operator": self.soc_operator,
            "scalar_relativistic_method": self.scalar_relativistic_method,
            "derivative_method": self.derivative_method,
            "active_space": self.active_space,
            "molecule_name": self.molecule_name,
            "atom_symbols": list(self.atom_symbols),
            "isotope_masses_amu": [float(value) for value in self.isotope_masses_amu],
            "reference_geometry_bohr": [
                [float(value) for value in row]
                for row in self.reference_geometry_bohr
            ],
            "calculation_input_sha256": self.calculation_input_sha256,
            "environment_sha256": self.environment_sha256,
            "traceable_nuclear_identity": self.traceable_nuclear_identity,
            "geometry_unit": self.geometry_unit,
            "energy_unit": self.energy_unit,
            "derivative_unit": self.derivative_unit,
            "extra": _canonical_v230(self.extra),
        }


def _validate_ladder_v230(name, levels, changes, tolerance):
    levels = tuple(str(value).strip() for value in levels)
    changes = tuple(float(value) for value in changes)
    tolerance = _optional_nonnegative_v230(f"{name}_tolerance", tolerance)
    supplied = bool(levels or changes or tolerance is not None)
    if not supplied:
        return levels, changes, tolerance, False
    if len(levels) < 2 or any(not value for value in levels):
        raise ValueError(f"{name}_levels must contain at least two nonempty labels.")
    if len(changes) != len(levels) - 1:
        raise ValueError(f"{name}_changes must contain one value between each level.")
    if not np.all(np.isfinite(changes)) or any(value < 0.0 for value in changes):
        raise ValueError(f"{name}_changes must be finite and nonnegative.")
    if tolerance is None:
        raise ValueError(f"{name}_tolerance is required with a convergence ladder.")
    return levels, changes, tolerance, bool(changes[-1] <= tolerance)


@dataclass(frozen=True)
class MolecularSOCValidationEvidenceV230:
    independent_reference_id: str | None = None
    independent_reference_error: float | None = None
    independent_reference_tolerance: float | None = None
    basis_levels: tuple[str, ...] = ()
    basis_changes: tuple[float, ...] = ()
    basis_tolerance: float | None = None
    method_levels: tuple[str, ...] = ()
    method_changes: tuple[float, ...] = ()
    method_tolerance: float | None = None
    translation_residual: float | None = None
    rotation_residual: float | None = None
    frame_invariance_tolerance: float | None = None
    tracking_minimum_overlap: float | None = None
    tracking_minimum_margin: float | None = None
    tracking_overlap_threshold: float | None = None
    tracking_margin_threshold: float | None = None

    def validate(self):
        reference_fields = (
            self.independent_reference_id,
            self.independent_reference_error,
            self.independent_reference_tolerance,
        )
        if any(value is not None for value in reference_fields):
            if any(value is None for value in reference_fields):
                raise ValueError("independent reference evidence is incomplete.")
            _required_text_v230("independent_reference_id", self.independent_reference_id)
            _optional_nonnegative_v230(
                "independent_reference_error", self.independent_reference_error
            )
            _optional_nonnegative_v230(
                "independent_reference_tolerance",
                self.independent_reference_tolerance,
            )
        _validate_ladder_v230(
            "basis", self.basis_levels, self.basis_changes, self.basis_tolerance
        )
        _validate_ladder_v230(
            "method", self.method_levels, self.method_changes, self.method_tolerance
        )
        invariance = (
            self.translation_residual,
            self.rotation_residual,
            self.frame_invariance_tolerance,
        )
        if any(value is not None for value in invariance):
            if any(value is None for value in invariance):
                raise ValueError("translation/rotation invariance evidence is incomplete.")
            for name, value in zip(
                (
                    "translation_residual",
                    "rotation_residual",
                    "frame_invariance_tolerance",
                ),
                invariance,
            ):
                _optional_nonnegative_v230(name, value)
        tracking = (
            self.tracking_minimum_overlap,
            self.tracking_minimum_margin,
            self.tracking_overlap_threshold,
            self.tracking_margin_threshold,
        )
        if any(value is not None for value in tracking):
            if any(value is None for value in tracking):
                raise ValueError("state-tracking evidence is incomplete.")
            for name, value in zip(
                (
                    "tracking_minimum_overlap",
                    "tracking_minimum_margin",
                    "tracking_overlap_threshold",
                    "tracking_margin_threshold",
                ),
                tracking,
            ):
                result = float(value)
                if not np.isfinite(result) or result < 0.0 or result > 1.0:
                    raise ValueError(f"{name} must lie in [0,1].")
        return self

    @property
    def independent_reference_validated(self):
        self.validate()
        return bool(
            self.independent_reference_id is not None
            and float(self.independent_reference_error)
            <= float(self.independent_reference_tolerance)
        )

    @property
    def basis_converged(self):
        return _validate_ladder_v230(
            "basis", self.basis_levels, self.basis_changes, self.basis_tolerance
        )[3]

    @property
    def method_converged(self):
        return _validate_ladder_v230(
            "method", self.method_levels, self.method_changes, self.method_tolerance
        )[3]

    @property
    def frame_invariance_validated(self):
        self.validate()
        return bool(
            self.translation_residual is not None
            and max(float(self.translation_residual), float(self.rotation_residual))
            <= float(self.frame_invariance_tolerance)
        )

    @property
    def state_tracking_validated(self):
        self.validate()
        return bool(
            self.tracking_minimum_overlap is not None
            and float(self.tracking_minimum_overlap)
            >= float(self.tracking_overlap_threshold)
            and float(self.tracking_minimum_margin)
            >= float(self.tracking_margin_threshold)
        )

    def as_dict(self):
        self.validate()
        return {
            "independent_reference_id": self.independent_reference_id,
            "independent_reference_error": self.independent_reference_error,
            "independent_reference_tolerance": self.independent_reference_tolerance,
            "independent_reference_validated": self.independent_reference_validated,
            "basis_levels": list(self.basis_levels),
            "basis_changes": [float(value) for value in self.basis_changes],
            "basis_tolerance": self.basis_tolerance,
            "basis_converged": self.basis_converged,
            "method_levels": list(self.method_levels),
            "method_changes": [float(value) for value in self.method_changes],
            "method_tolerance": self.method_tolerance,
            "method_converged": self.method_converged,
            "translation_residual": self.translation_residual,
            "rotation_residual": self.rotation_residual,
            "frame_invariance_tolerance": self.frame_invariance_tolerance,
            "frame_invariance_validated": self.frame_invariance_validated,
            "tracking_minimum_overlap": self.tracking_minimum_overlap,
            "tracking_minimum_margin": self.tracking_minimum_margin,
            "tracking_overlap_threshold": self.tracking_overlap_threshold,
            "tracking_margin_threshold": self.tracking_margin_threshold,
            "state_tracking_validated": self.state_tracking_validated,
        }


@dataclass(frozen=True)
class MolecularSOCAdmissionContractV230:
    capabilities: MolecularSOCCapabilitiesV230
    identity: MolecularSOCBackendIdentityV230
    evidence: MolecularSOCValidationEvidenceV230 = MolecularSOCValidationEvidenceV230()
    state_tracking_policy: str = "overlap_tracked_complete_multiplets"
    coordinate_definition: str = "explicit generalized coordinates in bohr"
    all_electronic_calculations_converged: bool = False

    def validate(self, symmetry_contract=None):
        self.capabilities.validate()
        self.identity.validate()
        self.evidence.validate()
        _required_text_v230("state_tracking_policy", self.state_tracking_policy)
        _required_text_v230("coordinate_definition", self.coordinate_definition)
        _native_bool_v230(
            "all_electronic_calculations_converged",
            self.all_electronic_calculations_converged,
        )
        if symmetry_contract is not None:
            if not isinstance(symmetry_contract, SOCSymmetryContractV221):
                raise TypeError("molecular admission requires SOCSymmetryContractV221.")
            if self.identity.electron_parity != symmetry_contract.electron_parity:
                raise ValueError(
                    "backend electron count disagrees with SOC electron parity."
                )
        return self

    @property
    def real_backend_admission_ready(self):
        self.validate()
        return bool(
            self.identity.real_ab_initio_source
            and self.identity.traceable_nuclear_identity
            and self.capabilities.trajectory_ready
            and self.all_electronic_calculations_converged
            and self.evidence.independent_reference_validated
            and self.evidence.basis_converged
            and self.evidence.method_converged
            and self.evidence.frame_invariance_validated
            and self.evidence.state_tracking_validated
        )

    def as_dict(self):
        self.validate()
        return {
            "capabilities": self.capabilities.as_dict(),
            "identity": self.identity.as_dict(),
            "evidence": self.evidence.as_dict(),
            "state_tracking_policy": self.state_tracking_policy,
            "coordinate_definition": self.coordinate_definition,
            "all_electronic_calculations_converged": bool(
                self.all_electronic_calculations_converged
            ),
            "real_backend_admission_ready": self.real_backend_admission_ready,
        }

    def fingerprint(self):
        payload = json.dumps(
            _canonical_v230(self.as_dict()),
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()


def provenance_with_molecular_soc_contract_v230(base_provenance, contract):
    """Bind the complete molecular admission declaration into operator identity."""
    base = base_provenance.validate()
    contract = contract.validate()
    parameters = dict(base.parameters)
    parameters.update(
        {
            "v230_molecular_soc_contract": contract.as_dict(),
            "v230_molecular_soc_contract_fingerprint": contract.fingerprint(),
        }
    )
    return ElectronicOperatorProvenanceV213(
        model_name=f"{base.model_name} via {contract.identity.backend_name}",
        model_version="v0.23.0-1",
        model_space=base.model_space,
        spin_free_method=contract.identity.electronic_method,
        soc_enabled=base.soc_enabled,
        soc_method=contract.identity.soc_operator,
        scalar_relativistic_method=contract.identity.scalar_relativistic_method,
        derivative_method=contract.identity.derivative_method,
        parameters=parameters,
    ).validate()


def molecular_soc_contract_from_provider_v230(provider):
    contract = getattr(provider, "molecular_soc_contract", None)
    if callable(contract):
        contract = contract()
    if not isinstance(contract, MolecularSOCAdmissionContractV230):
        raise TypeError(
            "provider must expose MolecularSOCAdmissionContractV230 as "
            "molecular_soc_contract."
        )
    symmetry = getattr(provider, "soc_symmetry_contract", None)
    if callable(symmetry):
        symmetry = symmetry()
    return contract.validate(symmetry)


def require_trajectory_ready_molecular_soc_v230(provider):
    contract = molecular_soc_contract_from_provider_v230(provider)
    if not contract.capabilities.trajectory_ready:
        raise ValueError(
            "moving-nuclear dynamics requires SOC derivatives, spin-free derivatives, "
            "derivative connections, and cross-geometry overlaps."
        )
    return contract
