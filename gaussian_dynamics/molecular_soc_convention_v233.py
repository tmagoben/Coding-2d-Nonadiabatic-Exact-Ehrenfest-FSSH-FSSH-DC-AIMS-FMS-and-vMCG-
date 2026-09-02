"""Frozen molecular spin-orbit matrix and derivative convention for v0.23.3."""

from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np

from .analytic_soc_models_v220 import SOCOperatorComponentsV220


SOC_CONVENTION_SCHEMA_V233 = "gnd-molecular-soc-matrix-convention-v0.23.3"
SOC_DERIVATIVE_SEMANTICS_V233 = (
    "physical fixed-frame operator derivative; moving-frame connection excluded"
)


def _canonical_bytes_v233(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class MolecularSOCMatrixConventionV233:
    schema: str
    operator_family: str
    one_electron_treatment: str
    two_electron_treatment: str
    mean_field_approximation: str
    prefactor_convention: str
    scalar_relativistic_method: str
    source_basis: str
    target_basis: str
    state_order: tuple
    electron_parity: str
    cartesian_component_order: tuple = ("x", "y", "z")
    spin_quantization_axis: str = "z"
    energy_unit: str = "hartree"
    coordinate_unit: str = "bohr"
    derivative_unit: str = "hartree/bohr"
    derivative_semantics: str = SOC_DERIVATIVE_SEMANTICS_V233
    complete_multiplets: bool = True
    external_magnetic_field: bool = False

    def validate(self):
        if self.schema != SOC_CONVENTION_SCHEMA_V233:
            raise ValueError("molecular SOC convention schema mismatch.")
        for name in (
            "operator_family",
            "one_electron_treatment",
            "two_electron_treatment",
            "mean_field_approximation",
            "prefactor_convention",
            "scalar_relativistic_method",
            "source_basis",
            "target_basis",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"SOC convention {name} cannot be empty.")
        if self.electron_parity not in {"even", "odd"}:
            raise ValueError("SOC convention electron parity must be even or odd.")
        if tuple(self.cartesian_component_order) != ("x", "y", "z"):
            raise ValueError("SOC Cartesian component order must be (x,y,z).")
        if self.spin_quantization_axis != "z":
            raise ValueError("v0.23.3 freezes the spin quantization axis to z.")
        if self.energy_unit != "hartree" or self.coordinate_unit != "bohr":
            raise ValueError("SOC convention requires hartree/bohr internal units.")
        if self.derivative_unit != "hartree/bohr":
            raise ValueError("SOC derivative unit must be hartree/bohr.")
        if self.derivative_semantics != SOC_DERIVATIVE_SEMANTICS_V233:
            raise ValueError("SOC derivative semantics are not the frozen physical form.")
        if type(self.complete_multiplets) is not bool or not self.complete_multiplets:
            raise ValueError("SOC convention requires complete multiplets.")
        if type(self.external_magnetic_field) is not bool:
            raise TypeError("external magnetic field flag must be a native Boolean.")
        if self.external_magnetic_field:
            raise ValueError("v0.23.3 molecular SOC convention excludes magnetic fields.")
        order = tuple(str(item) for item in self.state_order)
        if not order or len(set(order)) != len(order) or any(not item for item in order):
            raise ValueError("SOC convention state order must be unique and nonempty.")
        return self

    def as_dict(self):
        payload = asdict(self)
        payload["state_order"] = list(self.state_order)
        payload["cartesian_component_order"] = list(
            self.cartesian_component_order
        )
        return payload

    def fingerprint(self):
        self.validate()
        return hashlib.sha256(_canonical_bytes_v233(self.as_dict())).hexdigest()


def molecular_soc_convention_from_dict_v233(payload):
    if not isinstance(payload, dict):
        raise TypeError("molecular SOC convention payload must be a mapping.")
    expected = set(MolecularSOCMatrixConventionV233.__dataclass_fields__)
    if set(payload) != expected:
        raise ValueError("molecular SOC convention field set mismatch.")
    converted = dict(payload)
    converted["state_order"] = tuple(converted["state_order"])
    converted["cartesian_component_order"] = tuple(
        converted["cartesian_component_order"]
    )
    return MolecularSOCMatrixConventionV233(**converted).validate()


def analytic_soc_convention_v233(provider):
    provenance = provider.provenance.validate()
    parity = provider.soc_symmetry_contract.electron_parity
    return MolecularSOCMatrixConventionV233(
        schema=SOC_CONVENTION_SCHEMA_V233,
        operator_family=provenance.soc_method,
        one_electron_treatment="closed-form analytic model operator",
        two_electron_treatment="not applicable to analytic fixture",
        mean_field_approximation="none",
        prefactor_convention="all prefactors included in analytic model parameters",
        scalar_relativistic_method=provenance.scalar_relativistic_method,
        source_basis=provenance.model_space.representation,
        target_basis=provenance.model_space.representation,
        state_order=tuple(state.label for state in provenance.model_space.states),
        electron_parity=parity,
    ).validate()


@dataclass(frozen=True)
class MolecularSOCConventionAuditV233:
    convention_fingerprint: str
    component_hermiticity_residual: float
    derivative_hermiticity_residual: float
    checks: dict
    passed: bool

    def as_dict(self):
        return asdict(self)


def audit_molecular_soc_convention_v233(
    components,
    provenance,
    symmetry_contract,
    convention,
    *,
    tolerance=1.0e-12,
):
    if not isinstance(components, SOCOperatorComponentsV220):
        raise TypeError("SOC convention audit requires SOCOperatorComponentsV220.")
    components = components.validate()
    provenance = provenance.validate()
    convention = convention.validate()
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("SOC convention tolerance must be finite and positive.")
    state_labels = tuple(state.label for state in provenance.model_space.states)
    multiplicities = tuple(
        state.multiplicity for state in provenance.model_space.states
    )
    expected_parity = 0 if convention.electron_parity == "odd" else 1
    multiplicity_parity = bool(
        multiplicities
        and all(
            value is not None and int(value) % 2 == expected_parity
            for value in multiplicities
        )
    )
    H_residual = float(
        np.linalg.norm(components.H_soc - components.H_soc.conj().T, ord="fro")
    )
    K_residual = float(
        max(
            (
                np.linalg.norm(matrix - matrix.conj().T, ord="fro")
                for matrix in components.K_soc
            ),
            default=0.0,
        )
    )
    checks = {
        "exact_state_order": state_labels == tuple(convention.state_order),
        "complete_model_space": bool(provenance.model_space.complete_multiplets),
        "electron_parity": bool(
            symmetry_contract.electron_parity == convention.electron_parity
            and multiplicity_parity
        ),
        "zero_external_magnetic_field": not bool(
            symmetry_contract.external_magnetic_field
        ),
        "exact_scalar_relativistic_method": (
            provenance.scalar_relativistic_method
            == convention.scalar_relativistic_method
        ),
        "exact_soc_method": provenance.soc_method == convention.operator_family,
        "energy_unit": provenance.model_space.energy_unit == convention.energy_unit,
        "coordinate_unit": (
            provenance.model_space.coordinate_unit == convention.coordinate_unit
        ),
        "soc_hermiticity": H_residual <= tolerance,
        "soc_derivative_hermiticity": K_residual <= tolerance,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return MolecularSOCConventionAuditV233(
        convention_fingerprint=convention.fingerprint(),
        component_hermiticity_residual=H_residual,
        derivative_hermiticity_residual=K_residual,
        checks=checks,
        passed=bool(all(checks.values())),
    )


def require_molecular_soc_convention_v233(*args, **kwargs):
    report = audit_molecular_soc_convention_v233(*args, **kwargs)
    if not report.passed:
        failed = ", ".join(
            name for name, value in report.checks.items() if not value
        )
        raise ValueError("molecular SOC convention failed: " + failed)
    return report


def require_exact_molecular_soc_convention_v233(observed, trusted):
    """Require literal equality with a caller-supplied convention trust anchor."""
    if type(trusted) is not MolecularSOCMatrixConventionV233:
        raise TypeError("trusted convention must be MolecularSOCMatrixConventionV233.")
    if type(observed) is not MolecularSOCMatrixConventionV233:
        raise TypeError("observed convention must be MolecularSOCMatrixConventionV233.")
    trusted.validate()
    observed.validate()
    if observed != trusted:
        raise ValueError(
            "observed molecular SOC convention differs from the trusted "
            "operator, prefactor, basis, state-order, unit, or derivative identity."
        )
    return observed
