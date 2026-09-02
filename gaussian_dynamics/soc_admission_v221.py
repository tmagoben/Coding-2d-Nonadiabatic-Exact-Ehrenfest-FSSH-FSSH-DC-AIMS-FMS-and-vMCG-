"""Fail-closed SOC symmetry and provenance admission for v0.22.1."""

from dataclasses import asdict, dataclass
import numpy as np


def _scaled_error_v221(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        return float("inf")
    absolute = float(np.linalg.norm(left - right, ord="fro"))
    scale = max(
        float(np.linalg.norm(left, ord="fro")),
        float(np.linalg.norm(right, ord="fro")),
        1.0,
    )
    return absolute / scale


def _projector_residual_v221(projectors, nstate):
    if not isinstance(projectors, dict) or not projectors:
        return float("inf")
    matrices = []
    residuals = []
    for name, raw in projectors.items():
        if not isinstance(name, str) or not name.strip():
            return float("inf")
        projector = np.asarray(raw, dtype=complex)
        if projector.shape != (nstate, nstate) or not np.all(np.isfinite(projector)):
            return float("inf")
        matrices.append(projector)
        residuals.extend(
            (
                _scaled_error_v221(projector, projector.conj().T),
                _scaled_error_v221(projector @ projector, projector),
            )
        )
    for index, projector in enumerate(matrices):
        for other in matrices[index + 1 :]:
            residuals.append(float(np.linalg.norm(projector @ other, ord="fro")))
    residuals.append(
        _scaled_error_v221(
            sum(matrices, np.zeros((nstate, nstate), dtype=complex)),
            np.eye(nstate, dtype=complex),
        )
    )
    return max(residuals, default=0.0)


@dataclass(frozen=True)
class SOCSymmetryContractV221:
    electron_parity: str
    time_reversal_matrix: np.ndarray
    projectors: dict
    external_magnetic_field: bool = False

    def __post_init__(self):
        object.__setattr__(self, "electron_parity", str(self.electron_parity))
        object.__setattr__(
            self,
            "time_reversal_matrix",
            np.asarray(self.time_reversal_matrix, dtype=complex).copy(),
        )
        object.__setattr__(
            self,
            "projectors",
            {
                str(name): np.asarray(projector, dtype=complex).copy()
                for name, projector in dict(self.projectors).items()
            },
        )

    @property
    def fermionic(self):
        return self.electron_parity == "odd"

    def as_provenance_parameters(self):
        return {
            "electron_parity": self.electron_parity,
            "external_magnetic_field": bool(self.external_magnetic_field),
            "time_reversal_matrix": self.time_reversal_matrix.tolist(),
            "physical_projectors": {
                name: projector.tolist()
                for name, projector in sorted(self.projectors.items())
            },
        }


@dataclass(frozen=True)
class SOCSymmetryAuditV221:
    electron_parity: str
    nstate: int
    common_charge: int | None
    time_reversal_unitarity_residual: float
    time_reversal_square_residual: float
    projector_residual: float
    provenance_symmetry_residual: float
    checks: dict
    passed: bool
    tolerance: float

    def as_dict(self):
        return asdict(self)


def audit_soc_symmetry_contract_v221(
    model_space,
    contract,
    *,
    provenance=None,
    fermionic=None,
    tolerance=1.0e-12,
):
    model_space = model_space.validate()
    if not isinstance(contract, SOCSymmetryContractV221):
        raise TypeError("SOC symmetry admission requires SOCSymmetryContractV221.")
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("SOC symmetry tolerance must be finite and positive.")
    parity = contract.electron_parity
    if parity not in {"even", "odd"}:
        raise ValueError("electron_parity must be 'even' or 'odd'.")
    if not isinstance(contract.external_magnetic_field, (bool, np.bool_)):
        raise ValueError("external_magnetic_field must be Boolean.")
    expected_fermionic = parity == "odd"
    supplied_fermionic_matches = (
        True if fermionic is None else bool(fermionic) == expected_fermionic
    )

    multiplicities = [state.multiplicity for state in model_space.states]
    multiplicities_present = all(value is not None for value in multiplicities)
    expected_remainder = 0 if expected_fermionic else 1
    multiplicity_parity = multiplicities_present and all(
        int(value) % 2 == expected_remainder for value in multiplicities
    )
    charges = [state.charge for state in model_space.states]
    common_charge = (
        int(charges[0])
        if charges and charges[0] is not None and len(set(charges)) == 1
        else None
    )

    nstate = model_space.nstate
    J = np.asarray(contract.time_reversal_matrix, dtype=complex)
    if J.shape != (nstate, nstate) or not np.all(np.isfinite(J)):
        unitary_residual = float("inf")
        square_residual = float("inf")
    else:
        unitary_residual = _scaled_error_v221(
            J.conj().T @ J, np.eye(nstate, dtype=complex)
        )
        target = (-1.0 if expected_fermionic else 1.0) * np.eye(
            nstate, dtype=complex
        )
        square_residual = _scaled_error_v221(J @ J.conj(), target)
    projector_residual = _projector_residual_v221(contract.projectors, nstate)

    provenance_residuals = []
    provenance_keys_present = True
    if provenance is not None:
        provenance = provenance.validate()
        parameters = provenance.parameters
        required = {
            "electron_parity",
            "external_magnetic_field",
            "time_reversal_matrix",
            "physical_projectors",
        }
        provenance_keys_present = required <= set(parameters)
        if provenance_keys_present:
            provenance_residuals.append(
                0.0 if parameters["electron_parity"] == parity else float("inf")
            )
            provenance_residuals.append(
                0.0
                if bool(parameters["external_magnetic_field"])
                == bool(contract.external_magnetic_field)
                else float("inf")
            )
            provenance_residuals.append(
                _scaled_error_v221(parameters["time_reversal_matrix"], J)
            )
            stored_projectors = parameters["physical_projectors"]
            if set(stored_projectors) != set(contract.projectors):
                provenance_residuals.append(float("inf"))
            else:
                provenance_residuals.extend(
                    _scaled_error_v221(
                        stored_projectors[name], contract.projectors[name]
                    )
                    for name in contract.projectors
                )
        else:
            provenance_residuals.append(float("inf"))
    provenance_residual = max(provenance_residuals, default=0.0)

    checks = {
        "complete_multiplets": bool(model_space.complete_multiplets),
        "single_electron_parity": bool(multiplicity_parity),
        "single_charge_sector": common_charge is not None,
        "fermionic_declaration": bool(supplied_fermionic_matches),
        "zero_external_magnetic_field": not bool(contract.external_magnetic_field),
        "time_reversal_unitarity": unitary_residual <= tolerance,
        "time_reversal_square": square_residual <= tolerance,
        "physical_projectors": projector_residual <= tolerance,
        "symmetry_provenance_identity": (
            provenance_keys_present and provenance_residual <= tolerance
            if provenance is not None
            else True
        ),
    }
    return SOCSymmetryAuditV221(
        electron_parity=parity,
        nstate=nstate,
        common_charge=common_charge,
        time_reversal_unitarity_residual=float(unitary_residual),
        time_reversal_square_residual=float(square_residual),
        projector_residual=float(projector_residual),
        provenance_symmetry_residual=float(provenance_residual),
        checks={name: bool(value) for name, value in checks.items()},
        passed=bool(all(checks.values())),
        tolerance=tolerance,
    )


def soc_symmetry_contract_from_provider_v221(provider):
    if hasattr(provider, "soc_symmetry_contract"):
        contract = provider.soc_symmetry_contract
        if callable(contract):
            contract = contract()
    else:
        parameters = provider.provenance.parameters
        contract = SOCSymmetryContractV221(
            electron_parity=parameters.get("electron_parity", ""),
            time_reversal_matrix=provider.time_reversal_matrix,
            projectors=provider.projectors,
            external_magnetic_field=parameters.get(
                "external_magnetic_field", False
            ),
        )
    if not isinstance(contract, SOCSymmetryContractV221):
        raise TypeError("provider returned an invalid SOC symmetry contract.")
    return contract


def require_soc_symmetry_contract_v221(*args, **kwargs):
    report = audit_soc_symmetry_contract_v221(*args, **kwargs)
    if not report.passed:
        failed = ", ".join(
            name for name, passed in report.checks.items() if not passed
        )
        raise ValueError(f"SOC symmetry contract failed: {failed}.")
    return report
