"""Direct PySCF BP-SOMF state-interaction spin-orbit matrices for v0.24.1.

This module implements the static molecular-SOC tier only.  It deliberately does
not manufacture nuclear derivatives, derivative connections, or cross-geometry
overlaps.  Consequently, :class:`PySCFStateInteractionSOCProviderV241` cannot be
used as a moving-nuclei provider and cannot satisfy trajectory-ready admission.

The implementation uses PySCF's public scalar-orbital integral and FCI APIs.  The
spin algebra and assembly below are an independent implementation of the
Wigner--Eckart state-interaction construction.  No third-party SOC implementation
is bundled or imported.
"""

from dataclasses import asdict, dataclass, field, replace
import hashlib
from importlib import metadata
import importlib
import importlib.util
import json
import math

import numpy as np

from .electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from .molecular_soc_contract_v230 import (
    MolecularSOCAdmissionContractV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCCapabilitiesV230,
)
from .molecular_soc_convention_v233 import (
    MolecularSOCMatrixConventionV233,
    SOC_CONVENTION_SCHEMA_V233,
)
from .pyscf_nac_convention_v232 import PYSCF_REQUIRED_VERSION_V232
from .soc_admission_v221 import SOCSymmetryContractV221


PYSCF_REQUIRED_VERSION_V241 = PYSCF_REQUIRED_VERSION_V232
PYSCF_BP_SOMF_PROVIDER_NAME_V241 = "PySCFStateInteractionSOCProviderV241"
PYSCF_BP_SOMF_PROVIDER_VERSION_V241 = "0.24.1"
BP_SOMF_OPERATOR_FAMILY_V241 = "Breit-Pauli spin-orbit mean-field (BP-SOMF)"
BP_SOMF_ONE_ELECTRON_INTEGRAL_V241 = "int1e_prinvxp"
BP_SOMF_TWO_ELECTRON_INTEGRAL_V241 = "int2e_p1vxp1"
BP_SOMF_CARTESIAN_ORDER_V241 = ("x", "y", "z")
BP_SOMF_STATIC_LIMITATION_V241 = (
    "static SOC only: physical SOC derivatives, derivative connections, and "
    "cross-geometry many-electron overlaps are not implemented"
)


def _canonical_v241(value):
    if isinstance(value, np.generic):
        return _canonical_v241(value.item())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("v0.24.1 SOC data cannot contain non-finite values.")
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return _canonical_v241(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("v0.24.1 SOC dictionary keys must be strings.")
        return {
            key: _canonical_v241(item) for key, item in sorted(value.items())
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v241(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("v0.24.1 SOC data cannot contain non-finite values.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported v0.24.1 SOC value {type(value).__name__}.")


def _canonical_bytes_v241(value):
    return json.dumps(
        _canonical_v241(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _sha256_v241(value):
    return hashlib.sha256(_canonical_bytes_v241(value)).hexdigest()


def _finite_matrix_v241(name, value, *, ndim=None):
    matrix = np.asarray(value)
    if ndim is not None and matrix.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} contains non-finite data.")
    return matrix


def _complex_pairs_v241(value):
    array = np.asarray(value, dtype=complex)
    return {
        "real": np.asarray(array.real, dtype=float).tolist(),
        "imag": np.asarray(array.imag, dtype=float).tolist(),
    }


def _scaled_frobenius_error_v241(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        return float("inf")
    scale = max(
        float(np.linalg.norm(left, ord="fro")),
        float(np.linalg.norm(right, ord="fro")),
        1.0,
    )
    return float(np.linalg.norm(left - right, ord="fro") / scale)


@dataclass(frozen=True)
class PySCFStaticSOCProbeV241:
    installed: bool
    required_version: str
    distribution_version: str | None
    module_version: str | None
    exact_version: bool
    integral_apis_available: bool
    transition_rdm_api_available: bool
    spin_ladder_apis_available: bool
    usable: bool
    failure_reason: str | None

    def as_dict(self):
        return asdict(self)


def probe_pyscf_static_soc_runtime_v241():
    """Probe exactly the APIs used by the static BP-SOMF implementation."""

    try:
        installed = importlib.util.find_spec("pyscf") is not None
    except (ImportError, AttributeError, ValueError) as exc:
        return PySCFStaticSOCProbeV241(
            False,
            PYSCF_REQUIRED_VERSION_V241,
            None,
            None,
            False,
            False,
            False,
            False,
            False,
            f"PySCF discovery failed: {type(exc).__name__}: {exc}",
        )
    if not installed:
        return PySCFStaticSOCProbeV241(
            False,
            PYSCF_REQUIRED_VERSION_V241,
            None,
            None,
            False,
            False,
            False,
            False,
            False,
            "PySCF is not installed.",
        )

    try:
        distribution_version = metadata.version("pyscf")
    except metadata.PackageNotFoundError:
        distribution_version = None
    try:
        pyscf = importlib.import_module("pyscf")
        module_version = str(getattr(pyscf, "__version__", "unknown"))
        moleintor = importlib.import_module("pyscf.gto.moleintor")
        direct_spin1 = importlib.import_module("pyscf.fci.direct_spin1")
        addons = importlib.import_module("pyscf.fci.addons")
    except Exception as exc:
        return PySCFStaticSOCProbeV241(
            True,
            PYSCF_REQUIRED_VERSION_V241,
            distribution_version,
            None,
            False,
            False,
            False,
            False,
            False,
            f"PySCF SOC dependencies failed to import: {type(exc).__name__}: {exc}",
        )

    exact = bool(
        distribution_version == PYSCF_REQUIRED_VERSION_V241
        and module_version == PYSCF_REQUIRED_VERSION_V241
    )
    integral_table = getattr(moleintor, "_INTOR_FUNCTIONS", {})
    integral_apis = bool(
        BP_SOMF_ONE_ELECTRON_INTEGRAL_V241 in integral_table
        and BP_SOMF_TWO_ELECTRON_INTEGRAL_V241 in integral_table
    )
    transition_api = callable(getattr(direct_spin1, "trans_rdm1s", None))
    ladder_apis = all(
        callable(getattr(addons, name, None))
        for name in ("des_a", "des_b", "cre_a", "cre_b")
    )
    if not exact:
        reason = (
            f"PySCF {PYSCF_REQUIRED_VERSION_V241} is required; "
            f"distribution={distribution_version!r}, module={module_version!r}."
        )
    elif not integral_apis:
        reason = "PySCF lacks the required Breit-Pauli AO integral APIs."
    elif not transition_api:
        reason = "PySCF lacks spin-separated transition 1-RDM support."
    elif not ladder_apis:
        reason = "PySCF lacks the determinant spin-ladder primitives."
    else:
        reason = None
    usable = bool(exact and integral_apis and transition_api and ladder_apis)
    return PySCFStaticSOCProbeV241(
        True,
        PYSCF_REQUIRED_VERSION_V241,
        distribution_version,
        module_version,
        exact,
        integral_apis,
        transition_api,
        ladder_apis,
        usable,
        reason,
    )


def require_pyscf_static_soc_runtime_v241():
    probe = probe_pyscf_static_soc_runtime_v241()
    if not probe.installed:
        raise ImportError(
            "PySCF 2.13.1 is not installed; v0.24.1 static SOC fails closed."
        )
    if not probe.usable:
        raise RuntimeError(probe.failure_reason)
    return probe


def _twice_quantum_number_v241(name, value):
    if not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer twice-quantum-number.")
    return int(value)


def _factorial_integer_v241(value):
    rounded = int(round(value))
    if rounded < 0 or abs(float(value) - rounded) > 1.0e-12:
        raise ValueError("Clebsch-Gordan factorial argument is not nonnegative integer.")
    return math.factorial(rounded)


def clebsch_gordan_twice_v241(
    j1_twice,
    m1_twice,
    j2_twice,
    m2_twice,
    total_j_twice,
    total_m_twice,
):
    """Return ``<j1,m1;j2,m2|J,M>`` using integer twice-quantum numbers.

    The finite Racah sum avoids a SymPy runtime dependency and is adequate for the
    rank-one spin tensors used here.  Impossible angular-momentum combinations return
    exactly zero.
    """

    j1_twice = _twice_quantum_number_v241("2*j1", j1_twice)
    m1_twice = _twice_quantum_number_v241("2*m1", m1_twice)
    j2_twice = _twice_quantum_number_v241("2*j2", j2_twice)
    m2_twice = _twice_quantum_number_v241("2*m2", m2_twice)
    total_j_twice = _twice_quantum_number_v241("2*J", total_j_twice)
    total_m_twice = _twice_quantum_number_v241("2*M", total_m_twice)
    if min(j1_twice, j2_twice, total_j_twice) < 0:
        raise ValueError("angular momentum cannot be negative.")
    if m1_twice + m2_twice != total_m_twice:
        return 0.0
    angular_pairs = (
        (j1_twice, m1_twice),
        (j2_twice, m2_twice),
        (total_j_twice, total_m_twice),
    )
    if any(abs(m) > j or (j - m) % 2 for j, m in angular_pairs):
        return 0.0
    if (
        total_j_twice > j1_twice + j2_twice
        or total_j_twice < abs(j1_twice - j2_twice)
        or (j1_twice + j2_twice + total_j_twice) % 2
    ):
        return 0.0

    j1, m1, j2, m2, total_j, total_m = (
        value / 2.0
        for value in (
            j1_twice,
            m1_twice,
            j2_twice,
            m2_twice,
            total_j_twice,
            total_m_twice,
        )
    )
    triangle = (
        _factorial_integer_v241(total_j + j1 - j2)
        * _factorial_integer_v241(total_j - j1 + j2)
        * _factorial_integer_v241(j1 + j2 - total_j)
        / _factorial_integer_v241(j1 + j2 + total_j + 1.0)
    )
    projection = (
        _factorial_integer_v241(total_j + total_m)
        * _factorial_integer_v241(total_j - total_m)
        * _factorial_integer_v241(j1 - m1)
        * _factorial_integer_v241(j1 + m1)
        * _factorial_integer_v241(j2 - m2)
        * _factorial_integer_v241(j2 + m2)
    )
    prefactor = math.sqrt((2.0 * total_j + 1.0) * triangle * projection)
    lower = max(
        0,
        int(round(j2 - total_j - m1)),
        int(round(j1 + m2 - total_j)),
    )
    upper = min(
        int(round(j1 + j2 - total_j)),
        int(round(j1 - m1)),
        int(round(j2 + m2)),
    )
    total = 0.0
    for index in range(lower, upper + 1):
        denominator = (
            _factorial_integer_v241(index)
            * _factorial_integer_v241(j1 + j2 - total_j - index)
            * _factorial_integer_v241(j1 - m1 - index)
            * _factorial_integer_v241(j2 + m2 - index)
            * _factorial_integer_v241(total_j - j2 + m1 + index)
            * _factorial_integer_v241(total_j - j1 - m2 + index)
        )
        total += (-1.0) ** index / denominator
    return float(prefactor * total)


@dataclass(frozen=True)
class SpinFreeRootV241:
    label: str
    energy_hartree: float
    spin_twice: int
    reference_ms_twice: int
    spin_square: float | None = None

    def validate(self):
        if not str(self.label).strip():
            raise ValueError("spin-free root label cannot be empty.")
        energy = float(self.energy_hartree)
        if not np.isfinite(energy):
            raise ValueError("spin-free root energy must be finite.")
        spin_twice = _twice_quantum_number_v241("2*S", self.spin_twice)
        ms_twice = _twice_quantum_number_v241("2*M_S", self.reference_ms_twice)
        if spin_twice < 0 or abs(ms_twice) > spin_twice:
            raise ValueError("reference M_S lies outside the root spin multiplet.")
        if (spin_twice - ms_twice) % 2:
            raise ValueError("S and M_S must share integer/half-integer parity.")
        if self.spin_square is not None:
            observed = float(self.spin_square)
            expected = 0.25 * spin_twice * (spin_twice + 2)
            if not np.isfinite(observed) or abs(observed - expected) > 1.0e-6:
                raise ValueError("root <S^2> disagrees with its declared spin.")
        return self

    @property
    def multiplicity(self):
        return int(self.spin_twice) + 1

    def as_dict(self):
        self.validate()
        return {
            "label": self.label,
            "energy_hartree": float(self.energy_hartree),
            "spin_twice": int(self.spin_twice),
            "spin": 0.5 * int(self.spin_twice),
            "multiplicity": self.multiplicity,
            "reference_ms_twice": int(self.reference_ms_twice),
            "reference_ms": 0.5 * int(self.reference_ms_twice),
            "spin_square": None
            if self.spin_square is None
            else float(self.spin_square),
        }


def _format_half_integer_v241(value_twice):
    value_twice = int(value_twice)
    if value_twice % 2 == 0:
        return f"{value_twice // 2:+d}"
    sign = "+" if value_twice > 0 else "-"
    return f"{sign}{abs(value_twice)}/2"


@dataclass(frozen=True)
class SpinMicrostateV241:
    root_index: int
    root_label: str
    spin_twice: int
    ms_twice: int

    def validate(self):
        if int(self.root_index) != self.root_index or int(self.root_index) < 0:
            raise ValueError("microstate root index must be a nonnegative integer.")
        if not str(self.root_label).strip():
            raise ValueError("microstate root label cannot be empty.")
        spin_twice = _twice_quantum_number_v241("2*S", self.spin_twice)
        ms_twice = _twice_quantum_number_v241("2*M_S", self.ms_twice)
        if spin_twice < 0 or abs(ms_twice) > spin_twice:
            raise ValueError("microstate M_S lies outside its spin multiplet.")
        if (spin_twice - ms_twice) % 2:
            raise ValueError("microstate S/M_S parity mismatch.")
        return self

    @property
    def label(self):
        self.validate()
        return f"{self.root_label}(M={_format_half_integer_v241(self.ms_twice)})"

    @property
    def multiplicity(self):
        return int(self.spin_twice) + 1

    def as_dict(self):
        return {
            "root_index": int(self.root_index),
            "root_label": self.root_label,
            "spin_twice": int(self.spin_twice),
            "spin": 0.5 * int(self.spin_twice),
            "ms_twice": int(self.ms_twice),
            "ms": 0.5 * int(self.ms_twice),
            "multiplicity": self.multiplicity,
            "label": self.label,
        }


def complete_spin_microstates_v241(roots):
    """Expand roots into complete multiplets, ordered by root then decreasing M_S."""

    roots = tuple(root.validate() for root in roots)
    if not roots:
        raise ValueError("at least one spin-free root is required.")
    if len({root.label for root in roots}) != len(roots):
        raise ValueError("spin-free root labels must be unique.")
    parity = {root.multiplicity % 2 for root in roots}
    if len(parity) != 1:
        raise ValueError(
            "one molecular calculation cannot mix even- and odd-electron multiplets."
        )
    states = []
    for root_index, root in enumerate(roots):
        for ms_twice in range(root.spin_twice, -root.spin_twice - 1, -2):
            states.append(
                SpinMicrostateV241(
                    root_index,
                    root.label,
                    root.spin_twice,
                    ms_twice,
                ).validate()
            )
    return tuple(states)


def time_reversal_matrix_v241(microstates):
    """Unitary part of time reversal for the declared microstate order."""

    microstates = tuple(state.validate() for state in microstates)
    index = {
        (state.root_index, state.spin_twice, state.ms_twice): position
        for position, state in enumerate(microstates)
    }
    if len(index) != len(microstates):
        raise ValueError("microstate list contains duplicates.")
    result = np.zeros((len(microstates), len(microstates)), dtype=complex)
    for source, state in enumerate(microstates):
        target_key = (state.root_index, state.spin_twice, -state.ms_twice)
        if target_key not in index:
            raise ValueError("time reversal requires complete spin multiplets.")
        phase_exponent = (state.spin_twice + state.ms_twice) // 2
        result[index[target_key], source] = (-1.0) ** phase_exponent
    return result


def root_projectors_v241(microstates):
    microstates = tuple(state.validate() for state in microstates)
    labels = []
    for state in microstates:
        if state.root_label not in labels:
            labels.append(state.root_label)
    projectors = {}
    for label in labels:
        diagonal = [1.0 if state.root_label == label else 0.0 for state in microstates]
        projectors[label] = np.diag(diagonal).astype(complex)
    return projectors


@dataclass(frozen=True)
class BPSOMFIntegralsV241:
    one_electron_ao_cartesian: np.ndarray
    two_electron_somf_ao_cartesian: np.ndarray
    effective_ao_cartesian: np.ndarray
    effective_mo_cartesian: np.ndarray
    state_average_density_mo: np.ndarray
    state_average_density_ao: np.ndarray
    light_speed_au: float
    prefactor: float
    cartesian_component_order: tuple = BP_SOMF_CARTESIAN_ORDER_V241
    one_electron_integral: str = BP_SOMF_ONE_ELECTRON_INTEGRAL_V241
    two_electron_integral: str = BP_SOMF_TWO_ELECTRON_INTEGRAL_V241

    def __post_init__(self):
        for name in (
            "one_electron_ao_cartesian",
            "two_electron_somf_ao_cartesian",
            "effective_ao_cartesian",
            "effective_mo_cartesian",
            "state_average_density_mo",
            "state_average_density_ao",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name)).copy())
        object.__setattr__(
            self, "cartesian_component_order", tuple(self.cartesian_component_order)
        )

    def validate(self, tolerance=1.0e-10):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("BP-SOMF integral tolerance must be finite and positive.")
        one = _finite_matrix_v241(
            "one-electron SOC integrals", self.one_electron_ao_cartesian, ndim=3
        )
        two = _finite_matrix_v241(
            "two-electron SOMF integrals", self.two_electron_somf_ao_cartesian, ndim=3
        )
        effective_ao = _finite_matrix_v241(
            "effective AO SOC integrals", self.effective_ao_cartesian, ndim=3
        )
        effective_mo = _finite_matrix_v241(
            "effective MO SOC integrals", self.effective_mo_cartesian, ndim=3
        )
        density_mo = _finite_matrix_v241(
            "state-average MO density", self.state_average_density_mo, ndim=2
        )
        density_ao = _finite_matrix_v241(
            "state-average AO density", self.state_average_density_ao, ndim=2
        )
        if tuple(self.cartesian_component_order) != BP_SOMF_CARTESIAN_ORDER_V241:
            raise ValueError("BP-SOMF Cartesian component order must be (x,y,z).")
        if self.one_electron_integral != BP_SOMF_ONE_ELECTRON_INTEGRAL_V241:
            raise ValueError("wrong one-electron Breit-Pauli integral identity.")
        if self.two_electron_integral != BP_SOMF_TWO_ELECTRON_INTEGRAL_V241:
            raise ValueError("wrong two-electron Breit-Pauli integral identity.")
        if one.shape != two.shape or one.shape != effective_ao.shape:
            raise ValueError("AO SOC component tensors have incompatible shapes.")
        if one.shape[0] != 3 or one.shape[1] != one.shape[2]:
            raise ValueError("AO SOC tensor must have shape (3,nao,nao).")
        if effective_mo.shape[0] != 3 or effective_mo.shape[1] != effective_mo.shape[2]:
            raise ValueError("MO SOC tensor must have shape (3,nmo,nmo).")
        if density_mo.shape != effective_mo.shape[1:]:
            raise ValueError("MO density and SOC orbital dimensions differ.")
        if density_ao.shape != one.shape[1:]:
            raise ValueError("AO density and SOC AO dimensions differ.")
        light_speed = float(self.light_speed_au)
        prefactor = float(self.prefactor)
        if not np.isfinite(light_speed) or light_speed <= 0.0:
            raise ValueError("light speed must be finite and positive.")
        if not np.isfinite(prefactor) or prefactor <= 0.0:
            raise ValueError("BP-SOMF prefactor must be finite and positive.")
        if abs(prefactor - 0.5 / light_speed**2) > 1.0e-16:
            raise ValueError("BP-SOMF prefactor is not exactly 0.5/c^2.")
        for name, tensor in (
            ("one-electron", one),
            ("two-electron SOMF", two),
            ("effective AO", effective_ao),
        ):
            residual = float(np.max(np.abs(tensor + tensor.swapaxes(-1, -2))))
            if residual > tolerance:
                raise ValueError(f"{name} SOC integrals are not antisymmetric.")
        for component in effective_mo:
            if _scaled_frobenius_error_v241(component, component.conj().T) > tolerance:
                raise ValueError("effective MO SOC integrals are not Hermitian.")
        if _scaled_frobenius_error_v241(density_mo, density_mo.conj().T) > tolerance:
            raise ValueError("state-average MO density is not Hermitian.")
        if _scaled_frobenius_error_v241(density_ao, density_ao.conj().T) > tolerance:
            raise ValueError("state-average AO density is not Hermitian.")
        return self

    @property
    def one_electron_antisymmetry_residual(self):
        return float(
            np.max(
                np.abs(
                    self.one_electron_ao_cartesian
                    + self.one_electron_ao_cartesian.swapaxes(-1, -2)
                )
            )
        )

    @property
    def two_electron_antisymmetry_residual(self):
        return float(
            np.max(
                np.abs(
                    self.two_electron_somf_ao_cartesian
                    + self.two_electron_somf_ao_cartesian.swapaxes(-1, -2)
                )
            )
        )

    def as_dict(self, *, include_matrices=False):
        self.validate()
        payload = {
            "operator_family": BP_SOMF_OPERATOR_FAMILY_V241,
            "one_electron_integral": self.one_electron_integral,
            "two_electron_integral": self.two_electron_integral,
            "cartesian_component_order": list(self.cartesian_component_order),
            "light_speed_au": float(self.light_speed_au),
            "prefactor": float(self.prefactor),
            "one_electron_antisymmetry_residual": (
                self.one_electron_antisymmetry_residual
            ),
            "two_electron_antisymmetry_residual": (
                self.two_electron_antisymmetry_residual
            ),
            "one_electron_ao_norm": float(
                np.linalg.norm(self.one_electron_ao_cartesian)
            ),
            "two_electron_somf_ao_norm": float(
                np.linalg.norm(self.two_electron_somf_ao_cartesian)
            ),
            "effective_mo_norm_hartree": float(
                np.linalg.norm(self.effective_mo_cartesian)
            ),
            "state_average_density_trace": float(
                np.trace(self.state_average_density_mo).real
            ),
        }
        if include_matrices:
            payload["one_electron_ao_cartesian"] = _complex_pairs_v241(
                self.one_electron_ao_cartesian
            )
            payload["two_electron_somf_ao_cartesian"] = _complex_pairs_v241(
                self.two_electron_somf_ao_cartesian
            )
            payload["effective_ao_cartesian"] = _complex_pairs_v241(
                self.effective_ao_cartesian
            )
            payload["effective_mo_cartesian"] = _complex_pairs_v241(
                self.effective_mo_cartesian
            )
            payload["state_average_density_mo"] = _complex_pairs_v241(
                self.state_average_density_mo
            )
            payload["state_average_density_ao"] = _complex_pairs_v241(
                self.state_average_density_ao
            )
        return payload


def build_pyscf_bp_somf_integrals_v241(
    mol,
    mo_coeff,
    state_average_density_mo,
    *,
    imaginary_tolerance=1.0e-12,
):
    """Build all-electron BP-SOMF integrals using PySCF 2.13.1 primitives."""

    require_pyscf_static_soc_runtime_v241()
    from pyscf.lib.parameters import LIGHT_SPEED

    coefficients = np.asarray(mo_coeff, dtype=complex)
    density_mo = np.asarray(state_average_density_mo, dtype=complex)
    if coefficients.ndim != 2 or coefficients.shape[0] != int(mol.nao_nr()):
        raise ValueError("MO coefficients have the wrong AO dimension.")
    if density_mo.shape != (coefficients.shape[1], coefficients.shape[1]):
        raise ValueError("state-average MO density has the wrong shape.")
    if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(density_mo)):
        raise ValueError("MO coefficients and density must be finite.")
    if np.max(np.abs(coefficients.imag)) > imaginary_tolerance:
        raise ValueError(
            "v0.24.1 BP-SOMF supports real scalar common orbitals only."
        )
    if np.max(np.abs(density_mo.imag)) > imaginary_tolerance:
        raise ValueError(
            "v0.24.1 BP-SOMF supports a real state-average scalar density only."
        )
    coefficients = coefficients.real
    density_mo = density_mo.real
    density_ao = coefficients @ density_mo @ coefficients.T

    nao = int(mol.nao_nr())
    one_electron = np.zeros((3, nao, nao), dtype=float)
    for atom_index in range(int(mol.natm)):
        with mol.with_rinv_origin(mol.atom_coord(atom_index)):
            one_electron += float(mol.atom_charge(atom_index)) * np.asarray(
                mol.intor(BP_SOMF_ONE_ELECTRON_INTEGRAL_V241, comp=3),
                dtype=float,
            )

    two_electron_tensor = np.asarray(
        mol.intor(BP_SOMF_TWO_ELECTRON_INTEGRAL_V241, comp=3),
        dtype=float,
    )
    if two_electron_tensor.shape != (3, nao, nao, nao, nao):
        raise ValueError("PySCF returned an unexpected two-electron SOC tensor.")
    coulomb_like = np.einsum(
        "xijkl,kl->xij", two_electron_tensor, density_ao, optimize=True
    )
    exchange_left = np.einsum(
        "xijkl,jk->xil", two_electron_tensor, density_ao, optimize=True
    )
    exchange_right = np.einsum(
        "xijkl,li->xkj", two_electron_tensor, density_ao, optimize=True
    )
    two_electron_somf = (
        coulomb_like - 1.5 * exchange_left - 1.5 * exchange_right
    )
    light_speed = float(LIGHT_SPEED)
    prefactor = 0.5 / light_speed**2
    effective_ao = prefactor * (one_electron - two_electron_somf)
    transformed_real = np.einsum(
        "pi,xpq,qj->xij",
        coefficients,
        effective_ao,
        coefficients,
        optimize=True,
    )
    effective_mo = -1j * transformed_real.astype(complex)
    return BPSOMFIntegralsV241(
        one_electron,
        two_electron_somf,
        effective_ao,
        effective_mo,
        density_mo,
        density_ao,
        light_speed,
        prefactor,
    ).validate()


def state_average_density_mo_from_pyscf_ci_v241(
    ci_vectors,
    *,
    ncore,
    ncas,
    nmo,
    nelecas,
    weights,
):
    require_pyscf_static_soc_runtime_v241()
    from pyscf.fci import direct_spin1

    ci_vectors = tuple(np.asarray(ci) for ci in ci_vectors)
    weights = np.asarray(weights, dtype=float)
    ncore, ncas, nmo = int(ncore), int(ncas), int(nmo)
    nelecas = tuple(int(item) for item in nelecas)
    if not ci_vectors or weights.shape != (len(ci_vectors),):
        raise ValueError("one state-average weight is required per CI root.")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("state-average weights must be finite and nonnegative.")
    if abs(float(np.sum(weights)) - 1.0) > 1.0e-12:
        raise ValueError("state-average weights must sum to one.")
    if ncore < 0 or ncas < 1 or nmo < ncore + ncas:
        raise ValueError("invalid core/active/total orbital dimensions.")
    density = np.zeros((nmo, nmo), dtype=complex)
    density[:ncore, :ncore] = 2.0 * np.eye(ncore)
    active_slice = slice(ncore, ncore + ncas)
    for weight, ci in zip(weights, ci_vectors):
        density_alpha, density_beta = direct_spin1.make_rdm1s(
            ci, ncas, nelecas
        )
        density[active_slice, active_slice] += float(weight) * (
            np.asarray(density_alpha) + np.asarray(density_beta)
        )
    if _scaled_frobenius_error_v241(density, density.conj().T) > 1.0e-11:
        raise ValueError("state-average CI density is not Hermitian.")
    return density


def _apply_spin_ladder_once_v241(ci, ncas, nelecas, *, raising):
    from pyscf.fci import addons

    neleca, nelecb = (int(nelecas[0]), int(nelecas[1]))
    ci = np.asarray(ci)
    if raising:
        if nelecb == 0 or neleca == ncas:
            return np.zeros_like(ci), (neleca + 1, nelecb - 1)
        output = None
        for orbital in range(ncas):
            reduced = addons.des_b(ci, ncas, (neleca, nelecb), orbital)
            term = addons.cre_a(
                reduced, ncas, (neleca, nelecb - 1), orbital
            )
            output = term if output is None else output + term
        return output, (neleca + 1, nelecb - 1)
    if neleca == 0 or nelecb == ncas:
        return np.zeros_like(ci), (neleca - 1, nelecb + 1)
    output = None
    for orbital in range(ncas):
        reduced = addons.des_a(ci, ncas, (neleca, nelecb), orbital)
        term = addons.cre_b(reduced, ncas, (neleca - 1, nelecb), orbital)
        output = term if output is None else output + term
    return output, (neleca - 1, nelecb + 1)


def spin_ladder_ci_v241(
    ci,
    *,
    ncas,
    nelecas,
    spin_twice,
    target_ms_twice,
    norm_tolerance=1.0e-8,
):
    """Move a normalized spin-pure PySCF CI vector to another M_S component."""

    require_pyscf_static_soc_runtime_v241()
    ncas = int(ncas)
    current_nelec = tuple(int(item) for item in nelecas)
    current_ms_twice = current_nelec[0] - current_nelec[1]
    spin_twice = _twice_quantum_number_v241("2*S", spin_twice)
    target_ms_twice = _twice_quantum_number_v241("target 2*M_S", target_ms_twice)
    if abs(target_ms_twice) > spin_twice or (
        spin_twice - target_ms_twice
    ) % 2:
        raise ValueError("target M_S is not in the declared spin multiplet.")
    if (target_ms_twice - current_ms_twice) % 2:
        raise ValueError("spin ladder changes M_S in integer steps only.")
    raw_state = np.asarray(ci)
    if np.iscomplexobj(raw_state) and np.max(np.abs(raw_state.imag)) > 1.0e-12:
        raise ValueError(
            "v0.24.1 spin-ladder extraction supports real scalar CI vectors only."
        )
    # PySCF's direct-spin determinant ladder helpers allocate real work arrays.
    # Preserve that supported dtype instead of passing complex arrays with zero
    # imaginary part through the C-level transition-density routines.
    state = np.asarray(raw_state.real, dtype=float).copy()
    initial_norm = float(np.linalg.norm(state))
    if not np.isfinite(initial_norm) or abs(initial_norm - 1.0) > norm_tolerance:
        raise ValueError("spin ladder requires a normalized CI vector.")
    spin = 0.5 * spin_twice
    while current_ms_twice != target_ms_twice:
        current_ms = 0.5 * current_ms_twice
        raising = target_ms_twice > current_ms_twice
        if raising:
            ladder_norm = math.sqrt(
                spin * (spin + 1.0) - current_ms * (current_ms + 1.0)
            )
        else:
            ladder_norm = math.sqrt(
                spin * (spin + 1.0) - current_ms * (current_ms - 1.0)
            )
        if ladder_norm <= 0.0:
            raise ValueError("requested spin-ladder operation annihilates this state.")
        state, current_nelec = _apply_spin_ladder_once_v241(
            state, ncas, current_nelec, raising=raising
        )
        state = np.asarray(state, dtype=float) / ladder_norm
        current_ms_twice += 2 if raising else -2
        observed_norm = float(np.linalg.norm(state))
        if not np.isfinite(observed_norm) or abs(observed_norm - 1.0) > norm_tolerance:
            raise ValueError(
                "PySCF CI spin-ladder normalization disagrees with declared S/M_S."
            )
    return state, current_nelec


def _allowed_sample_ms_v241(left, right):
    minimum_spin = min(left.spin_twice, right.spin_twice)
    candidates = list(range(minimum_spin, -minimum_spin - 1, -2))
    candidates.sort(
        key=lambda ms: (
            abs(ms - left.reference_ms_twice)
            + abs(ms - right.reference_ms_twice),
            -ms,
        )
    )
    for sample_ms in candidates:
        coefficient = clebsch_gordan_twice_v241(
            right.spin_twice,
            sample_ms,
            2,
            0,
            left.spin_twice,
            sample_ms,
        )
        if abs(coefficient) > 1.0e-14:
            return sample_ms, coefficient
    return None, 0.0


def wigner_reduced_transition_density_from_pyscf_ci_v241(
    ci_vectors,
    roots,
    *,
    ncore,
    ncas,
    nmo,
    nelecas,
):
    """Extract rank-one reduced spin transition densities from PySCF CI roots."""

    require_pyscf_static_soc_runtime_v241()
    from pyscf.fci import direct_spin1

    roots = tuple(root.validate() for root in roots)
    ci_vectors = tuple(np.asarray(ci) for ci in ci_vectors)
    ncore, ncas, nmo = int(ncore), int(ncas), int(nmo)
    nelecas = tuple(int(item) for item in nelecas)
    if len(ci_vectors) != len(roots) or not roots:
        raise ValueError("one CI vector is required per spin-free root.")
    reference_ms_twice = nelecas[0] - nelecas[1]
    if any(root.reference_ms_twice != reference_ms_twice for root in roots):
        raise ValueError("root reference M_S disagrees with PySCF nelecas.")
    electron_parities = {(root.spin_twice + 1) % 2 for root in roots}
    if len(electron_parities) != 1:
        raise ValueError("PySCF roots cannot mix electron-number parity sectors.")
    reduced = np.zeros((len(roots), len(roots), nmo, nmo), dtype=complex)
    sample_ms_table = np.zeros((len(roots), len(roots)), dtype=int)
    coefficient_table = np.zeros((len(roots), len(roots)), dtype=float)
    active = slice(ncore, ncore + ncas)
    cache = {}

    def lifted(index, target_ms):
        key = (index, target_ms)
        if key not in cache:
            cache[key] = spin_ladder_ci_v241(
                ci_vectors[index],
                ncas=ncas,
                nelecas=nelecas,
                spin_twice=roots[index].spin_twice,
                target_ms_twice=target_ms,
            )
        return cache[key]

    for left_index, left in enumerate(roots):
        for right_index, right in enumerate(roots):
            if (
                abs(left.spin_twice - right.spin_twice) > 2
                or left.spin_twice + right.spin_twice < 2
            ):
                continue
            sample_ms, coefficient = _allowed_sample_ms_v241(left, right)
            if sample_ms is None:
                raise ValueError(
                    "could not find a nonzero rank-one Wigner reference component."
                )
            left_ci, left_nelec = lifted(left_index, sample_ms)
            right_ci, right_nelec = lifted(right_index, sample_ms)
            if left_nelec != right_nelec:
                raise ValueError("Wigner reference CI vectors occupy different sectors.")
            density_alpha, density_beta = direct_spin1.trans_rdm1s(
                left_ci, right_ci, ncas, left_nelec
            )
            spin_density = np.asarray(density_alpha) - np.asarray(density_beta)
            reduced[left_index, right_index, active, active] = (
                spin_density / (math.sqrt(2.0) * coefficient)
            )
            sample_ms_table[left_index, right_index] = sample_ms
            coefficient_table[left_index, right_index] = coefficient
    if not np.all(np.isfinite(reduced)):
        raise ValueError("Wigner-reduced transition density contains non-finite data.")
    return reduced, sample_ms_table, coefficient_table


def _one_body_transition_contraction_v241(operator, density):
    """Contract PySCF's integral/RDM index convention without reindexing.

    ``direct_spin1.trans_rdm1s`` documents its returned entry ``D[p,q]`` as
    ``<bra|a_q^+ a_p|ket>``.  The ``prinvxp``/``p1vxp1`` spin-orbit integral
    convention is aligned with those density indices, so the contraction is over
    identical ``p,q`` positions.  This orientation is frozen and tested against an
    independent PySCF JK contraction path in the v0.24.1 runtime evidence.
    """

    return np.einsum("pq,...pq->...", operator, density, optimize=True)


@dataclass(frozen=True)
class StateInteractionSOCMatricesV241:
    roots: tuple
    microstates: tuple
    H_spin_free: np.ndarray
    H_soc: np.ndarray
    H_total: np.ndarray
    soc_eigenvalues_hartree: np.ndarray
    soc_eigenvectors: np.ndarray
    time_reversal_matrix: np.ndarray
    root_projectors: dict
    hermiticity_residual: float
    time_reversal_residual: float
    time_reversal_square_residual: float
    maximum_kramers_pair_splitting_hartree: float | None

    def __post_init__(self):
        object.__setattr__(self, "roots", tuple(self.roots))
        object.__setattr__(self, "microstates", tuple(self.microstates))
        for name in (
            "H_spin_free",
            "H_soc",
            "H_total",
            "soc_eigenvalues_hartree",
            "soc_eigenvectors",
            "time_reversal_matrix",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name)).copy())
        object.__setattr__(
            self,
            "root_projectors",
            {
                str(name): np.asarray(matrix, dtype=complex).copy()
                for name, matrix in dict(self.root_projectors).items()
            },
        )

    def validate(self, tolerance=1.0e-10):
        roots = tuple(root.validate() for root in self.roots)
        expected_microstates = complete_spin_microstates_v241(roots)
        if tuple(self.microstates) != expected_microstates:
            raise ValueError("state-interaction microstate order is not complete/exact.")
        nstate = len(expected_microstates)
        H0 = _finite_matrix_v241("H_spin_free", self.H_spin_free, ndim=2)
        Hso = _finite_matrix_v241("H_soc", self.H_soc, ndim=2)
        H = _finite_matrix_v241("H_total", self.H_total, ndim=2)
        eigenvalues = _finite_matrix_v241(
            "SOC eigenvalues", self.soc_eigenvalues_hartree, ndim=1
        )
        eigenvectors = _finite_matrix_v241(
            "SOC eigenvectors", self.soc_eigenvectors, ndim=2
        )
        time_reversal = _finite_matrix_v241(
            "time-reversal matrix", self.time_reversal_matrix, ndim=2
        )
        expected_shape = (nstate, nstate)
        if any(
            item.shape != expected_shape
            for item in (H0, Hso, H, eigenvectors, time_reversal)
        ) or eigenvalues.shape != (nstate,):
            raise ValueError("state-interaction matrices have incompatible dimensions.")
        if _scaled_frobenius_error_v241(H, H0 + Hso) > tolerance:
            raise ValueError("H_total is not H_spin_free + H_soc.")
        for name, matrix in (("H_spin_free", H0), ("H_soc", Hso), ("H_total", H)):
            if _scaled_frobenius_error_v241(matrix, matrix.conj().T) > tolerance:
                raise ValueError(f"{name} is not Hermitian.")
        identity = np.eye(nstate, dtype=complex)
        if _scaled_frobenius_error_v241(
            eigenvectors.conj().T @ eigenvectors, identity
        ) > tolerance:
            raise ValueError("SOC eigenvectors are not unitary.")
        reconstructed = (
            eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T
        )
        if _scaled_frobenius_error_v241(reconstructed, H) > tolerance:
            raise ValueError("SOC eigensystem does not reconstruct H_total.")
        if _scaled_frobenius_error_v241(
            time_reversal.conj().T @ time_reversal, identity
        ) > tolerance:
            raise ValueError("time-reversal matrix is not unitary.")
        odd = roots[0].multiplicity % 2 == 0
        target_square = (-1.0 if odd else 1.0) * identity
        if _scaled_frobenius_error_v241(
            time_reversal @ time_reversal.conj(), target_square
        ) > tolerance:
            raise ValueError("time-reversal operator has the wrong square.")
        if _scaled_frobenius_error_v241(
            H, time_reversal @ H.conj() @ time_reversal.conj().T
        ) > tolerance:
            raise ValueError("state-interaction Hamiltonian violates time reversal.")
        expected_projectors = root_projectors_v241(expected_microstates)
        if set(self.root_projectors) != set(expected_projectors):
            raise ValueError("root projectors do not match the spin-free roots.")
        for name in expected_projectors:
            if _scaled_frobenius_error_v241(
                self.root_projectors[name], expected_projectors[name]
            ) > tolerance:
                raise ValueError("root projector differs from the microstate order.")
        if abs(float(self.hermiticity_residual)) > tolerance:
            raise ValueError("recorded state-interaction Hermiticity residual failed.")
        if abs(float(self.time_reversal_residual)) > tolerance:
            raise ValueError("recorded state-interaction time-reversal residual failed.")
        if abs(float(self.time_reversal_square_residual)) > tolerance:
            raise ValueError("recorded time-reversal-square residual failed.")
        if odd:
            if self.maximum_kramers_pair_splitting_hartree is None:
                raise ValueError("odd-electron SOC requires a Kramers splitting metric.")
            if float(self.maximum_kramers_pair_splitting_hartree) > tolerance:
                raise ValueError("Kramers pairs are split above tolerance.")
        elif self.maximum_kramers_pair_splitting_hartree is not None:
            raise ValueError("even-electron SOC must not report Kramers pairing.")
        return self

    @property
    def state_order(self):
        return tuple(state.label for state in self.microstates)

    def as_dict(self, *, include_matrices=True):
        self.validate()
        payload = {
            "roots": [root.as_dict() for root in self.roots],
            "microstates": [state.as_dict() for state in self.microstates],
            "state_order": list(self.state_order),
            "hermiticity_residual": float(self.hermiticity_residual),
            "time_reversal_residual": float(self.time_reversal_residual),
            "time_reversal_square_residual": float(
                self.time_reversal_square_residual
            ),
            "maximum_kramers_pair_splitting_hartree": (
                self.maximum_kramers_pair_splitting_hartree
            ),
            "soc_eigenvalues_hartree": np.asarray(
                self.soc_eigenvalues_hartree, dtype=float
            ).tolist(),
        }
        if include_matrices:
            payload.update(
                {
                    "H_spin_free": _complex_pairs_v241(self.H_spin_free),
                    "H_soc": _complex_pairs_v241(self.H_soc),
                    "H_total": _complex_pairs_v241(self.H_total),
                    "soc_eigenvectors": _complex_pairs_v241(
                        self.soc_eigenvectors
                    ),
                    "time_reversal_matrix": _complex_pairs_v241(
                        self.time_reversal_matrix
                    ),
                    "root_projectors": {
                        name: _complex_pairs_v241(matrix)
                        for name, matrix in sorted(self.root_projectors.items())
                    },
                }
            )
        return payload


def assemble_state_interaction_soc_v241(
    roots,
    wigner_reduced_transition_density,
    effective_mo_cartesian,
    *,
    tolerance=1.0e-10,
):
    """Assemble direct H_sf, H_SOC, and H_total in a complete |root,S,M_S> basis."""

    roots = tuple(root.validate() for root in roots)
    microstates = complete_spin_microstates_v241(roots)
    reduced = np.asarray(wigner_reduced_transition_density, dtype=complex)
    integrals = np.asarray(effective_mo_cartesian, dtype=complex)
    nroot = len(roots)
    if reduced.ndim != 4 or reduced.shape[:2] != (nroot, nroot):
        raise ValueError("reduced transition densities have the wrong root dimensions.")
    if integrals.shape != (3, reduced.shape[2], reduced.shape[3]):
        raise ValueError("SOC integrals and transition densities have incompatible shapes.")
    if reduced.shape[2] != reduced.shape[3] or not np.all(np.isfinite(reduced)):
        raise ValueError("reduced transition densities must be finite square matrices.")
    if not np.all(np.isfinite(integrals)):
        raise ValueError("SOC integrals contain non-finite data.")

    h_x, h_y, h_z = integrals
    reduced_amplitudes = {
        0: _one_body_transition_contraction_v241(h_z, reduced) / math.sqrt(2.0),
        1: -0.5
        * _one_body_transition_contraction_v241(h_x - 1j * h_y, reduced),
        -1: 0.5
        * _one_body_transition_contraction_v241(h_x + 1j * h_y, reduced),
    }
    nstate = len(microstates)
    H_soc = np.zeros((nstate, nstate), dtype=complex)
    for row, bra in enumerate(microstates):
        for column, ket in enumerate(microstates):
            delta_ms_twice = bra.ms_twice - ket.ms_twice
            if delta_ms_twice not in {-2, 0, 2}:
                continue
            spherical_component = delta_ms_twice // 2
            coefficient = clebsch_gordan_twice_v241(
                ket.spin_twice,
                ket.ms_twice,
                2,
                2 * spherical_component,
                bra.spin_twice,
                bra.ms_twice,
            )
            H_soc[row, column] = coefficient * reduced_amplitudes[
                spherical_component
            ][bra.root_index, ket.root_index]
    hermiticity_residual = float(np.linalg.norm(H_soc - H_soc.conj().T, ord="fro"))
    if hermiticity_residual > tolerance:
        raise ValueError(
            "direct state-interaction H_SOC is not Hermitian; root phases, transition "
            "density orientation, spin labels, or operator convention are inconsistent."
        )

    energies = np.asarray(
        [roots[state.root_index].energy_hartree for state in microstates],
        dtype=float,
    )
    H_spin_free = np.diag(energies).astype(complex)
    H_total = H_spin_free + H_soc
    eigenvalues, eigenvectors = np.linalg.eigh(H_total)
    time_reversal = time_reversal_matrix_v241(microstates)
    time_reversal_residual = _scaled_frobenius_error_v241(
        H_total, time_reversal @ H_total.conj() @ time_reversal.conj().T
    )
    target_square = (
        -np.eye(nstate, dtype=complex)
        if roots[0].multiplicity % 2 == 0
        else np.eye(nstate, dtype=complex)
    )
    time_reversal_square_residual = _scaled_frobenius_error_v241(
        time_reversal @ time_reversal.conj(), target_square
    )
    if roots[0].multiplicity % 2 == 0:
        if nstate % 2:
            raise ValueError("odd-electron complete multiplets require even dimension.")
        kramers_splitting = float(
            np.max(np.abs(eigenvalues[1::2] - eigenvalues[0::2]))
        )
    else:
        kramers_splitting = None
    return StateInteractionSOCMatricesV241(
        roots=roots,
        microstates=microstates,
        H_spin_free=H_spin_free,
        H_soc=H_soc,
        H_total=H_total,
        soc_eigenvalues_hartree=eigenvalues,
        soc_eigenvectors=eigenvectors,
        time_reversal_matrix=time_reversal,
        root_projectors=root_projectors_v241(microstates),
        hermiticity_residual=hermiticity_residual,
        time_reversal_residual=time_reversal_residual,
        time_reversal_square_residual=time_reversal_square_residual,
        maximum_kramers_pair_splitting_hartree=kramers_splitting,
    ).validate(tolerance=tolerance)


@dataclass(frozen=True)
class PySCFStateInteractionSOCResultV241:
    matrices: StateInteractionSOCMatricesV241
    integrals: BPSOMFIntegralsV241
    wigner_reduced_transition_density: np.ndarray
    wigner_sample_ms_twice: np.ndarray
    wigner_reference_coefficients: np.ndarray
    capabilities: MolecularSOCCapabilitiesV230
    identity: MolecularSOCBackendIdentityV230
    molecular_soc_contract: MolecularSOCAdmissionContractV230
    convention: MolecularSOCMatrixConventionV233
    provenance: ElectronicOperatorProvenanceV213
    scf_converged: bool
    casscf_converged: bool
    soc_assembled: bool
    runtime_probe: PySCFStaticSOCProbeV241
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in (
            "wigner_reduced_transition_density",
            "wigner_sample_ms_twice",
            "wigner_reference_coefficients",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name)).copy())
        object.__setattr__(self, "metadata", dict(self.metadata))

    def validate(self, tolerance=1.0e-10):
        self.matrices.validate(tolerance=tolerance)
        self.integrals.validate(tolerance=tolerance)
        self.capabilities.validate()
        self.identity.validate()
        self.molecular_soc_contract.validate()
        self.convention.validate()
        self.provenance.validate()
        if self.molecular_soc_contract.capabilities != self.capabilities:
            raise ValueError("result capabilities and molecular contract differ.")
        if self.molecular_soc_contract.identity != self.identity:
            raise ValueError("result identity and molecular contract differ.")
        if self.capabilities.tier != "static_soc" or self.capabilities.trajectory_ready:
            raise ValueError("v0.24.1 PySCF BP-SOMF result must remain static-only.")
        if self.molecular_soc_contract.real_backend_admission_ready:
            raise ValueError("static PySCF SOC must not pass trajectory admission.")
        if tuple(self.convention.state_order) != self.matrices.state_order:
            raise ValueError("SOC convention and direct matrix state order differ.")
        if self.convention.operator_family != BP_SOMF_OPERATOR_FAMILY_V241:
            raise ValueError("result uses an unexpected SOC operator family.")
        if self.provenance.soc_method != BP_SOMF_OPERATOR_FAMILY_V241:
            raise ValueError("operator provenance and SOC convention differ.")
        if self.provenance.model_space.nstate != len(self.matrices.microstates):
            raise ValueError("operator provenance has the wrong model-space dimension.")
        if self.runtime_probe.usable is not True:
            raise ValueError("result lacks a usable exact PySCF runtime probe.")
        for name, value in (
            ("scf_converged", self.scf_converged),
            ("casscf_converged", self.casscf_converged),
            ("soc_assembled", self.soc_assembled),
        ):
            if type(value) is not bool or not value:
                raise ValueError(f"{name} must be an affirmative native Boolean.")
        reduced = _finite_matrix_v241(
            "Wigner-reduced transition density",
            self.wigner_reduced_transition_density,
            ndim=4,
        )
        samples = np.asarray(self.wigner_sample_ms_twice)
        coefficients = _finite_matrix_v241(
            "Wigner reference coefficients",
            self.wigner_reference_coefficients,
            ndim=2,
        )
        nroot = len(self.matrices.roots)
        nmo = self.integrals.effective_mo_cartesian.shape[1]
        if reduced.shape != (nroot, nroot, nmo, nmo):
            raise ValueError("Wigner density dimensions differ from roots/orbitals.")
        if samples.shape != (nroot, nroot) or coefficients.shape != (nroot, nroot):
            raise ValueError("Wigner reference tables have the wrong shape.")
        if not np.issubdtype(samples.dtype, np.integer):
            raise ValueError("Wigner sample M_S table must contain integers.")
        if not isinstance(self.metadata, dict):
            raise TypeError("SOC result metadata must be a dictionary.")
        _canonical_v241(self.metadata)
        return self

    @property
    def static_soc_admitted(self):
        self.validate()
        return True

    @property
    def trajectory_ready(self):
        return False

    def fingerprint(self):
        return _sha256_v241(self.as_dict(include_large_arrays=True))

    def as_dict(self, *, include_large_arrays=False):
        self.validate()
        payload = {
            "schema": "gnd-pyscf-bp-somf-state-interaction-v0.24.1",
            "provider": PYSCF_BP_SOMF_PROVIDER_NAME_V241,
            "provider_version": PYSCF_BP_SOMF_PROVIDER_VERSION_V241,
            "static_limitation": BP_SOMF_STATIC_LIMITATION_V241,
            "runtime_probe": self.runtime_probe.as_dict(),
            "capabilities": self.capabilities.as_dict(),
            "identity": self.identity.as_dict(),
            "molecular_soc_contract": self.molecular_soc_contract.as_dict(),
            "convention": self.convention.as_dict(),
            "convention_fingerprint": self.convention.fingerprint(),
            "provenance": self.provenance.as_dict(),
            "provenance_fingerprint": self.provenance.fingerprint(),
            "scf_converged": self.scf_converged,
            "casscf_converged": self.casscf_converged,
            "soc_assembled": self.soc_assembled,
            "static_soc_admitted": True,
            "trajectory_ready": False,
            "matrices": self.matrices.as_dict(include_matrices=True),
            "integrals": self.integrals.as_dict(
                include_matrices=include_large_arrays
            ),
            "wigner_sample_ms_twice": np.asarray(
                self.wigner_sample_ms_twice, dtype=int
            ).tolist(),
            "wigner_reference_coefficients": np.asarray(
                self.wigner_reference_coefficients, dtype=float
            ).tolist(),
            "metadata": _canonical_v241(self.metadata),
        }
        if include_large_arrays:
            payload["wigner_reduced_transition_density"] = _complex_pairs_v241(
                self.wigner_reduced_transition_density
            )
        return payload


def _model_space_v241(microstates, charge):
    return ElectronicModelSpaceV213(
        name="v0.24.1 PySCF BP-SOMF complete spin-microstate space",
        representation="fixed_spin_diabatic",
        states=tuple(
            ElectronicStateDescriptorV213(
                label=state.label,
                source_root=state.root_label,
                multiplicity=state.multiplicity,
                component=f"M={_format_half_integer_v241(state.ms_twice)}",
                charge=int(charge),
            )
            for state in microstates
        ),
        complete_multiplets=True,
    ).validate()


def _normalised_weights_v241(weights, nroot):
    if weights is None:
        return np.full(nroot, 1.0 / nroot)
    values = np.asarray(weights, dtype=float)
    if values.shape != (nroot,) or not np.all(np.isfinite(values)):
        raise ValueError("one finite state-average weight is required per root.")
    if np.any(values < 0.0) or abs(float(np.sum(values)) - 1.0) > 1.0e-12:
        raise ValueError("state-average weights must be nonnegative and sum to one.")
    return values


def _root_spin_data_v241(ci_vectors, *, ncas, nelecas, declared=None):
    from pyscf.fci import spin_op

    if declared is not None:
        declared = tuple(int(value) for value in declared)
        if len(declared) != len(ci_vectors):
            raise ValueError("one declared spin is required per CI root.")
    output = []
    for index, ci in enumerate(ci_vectors):
        spin_square, multiplicity = spin_op.spin_square(
            ci, int(ncas), tuple(nelecas)
        )
        observed_spin_twice = int(round(float(multiplicity) - 1.0))
        if declared is not None and observed_spin_twice != declared[index]:
            raise ValueError("declared root spin disagrees with PySCF <S^2>.")
        spin_twice = observed_spin_twice if declared is None else declared[index]
        expected = 0.25 * spin_twice * (spin_twice + 2)
        if abs(float(spin_square) - expected) > 1.0e-6:
            raise ValueError("PySCF CI root is not spin-pure at the declared S.")
        output.append((spin_twice, float(spin_square)))
    return tuple(output)


class PySCFStateInteractionSOCProviderV241:
    """Static direct BP-SOMF state-interaction provider for a converged PySCF CASSCF.

    ``components`` and moving-geometry methods intentionally raise.  This prevents a
    valid static matrix from being mistaken for the full trajectory provider contract.
    """

    adapter_name = PYSCF_BP_SOMF_PROVIDER_NAME_V241
    adapter_version = PYSCF_BP_SOMF_PROVIDER_VERSION_V241

    def __init__(
        self,
        mc,
        *,
        environment_sha256,
        root_labels=None,
        root_spin_twice=None,
        weights=None,
        molecule_name="PySCF molecule",
        basis_label=None,
        isotope_masses_amu=None,
    ):
        self.runtime_probe = require_pyscf_static_soc_runtime_v241()
        self.mc = mc
        self.mol = getattr(mc, "mol", None)
        if self.mol is None:
            raise TypeError("PySCF provider requires a CASSCF object with mol.")
        if type(getattr(getattr(mc, "_scf", None), "converged", None)) not in {
            bool,
            np.bool_,
        } or not bool(mc._scf.converged):
            raise RuntimeError("PySCF SCF reference is not converged.")
        if type(getattr(mc, "converged", None)) not in {bool, np.bool_} or not bool(
            mc.converged
        ):
            raise RuntimeError("PySCF CASSCF calculation is not converged.")
        if not isinstance(environment_sha256, str) or len(environment_sha256) != 64:
            raise ValueError("environment_sha256 must be a lowercase SHA-256 digest.")
        if any(character not in "0123456789abcdef" for character in environment_sha256):
            raise ValueError("environment_sha256 must be a lowercase SHA-256 digest.")

        raw_ci = getattr(mc, "ci", None)
        ci_vectors = tuple(raw_ci) if isinstance(raw_ci, (tuple, list)) else (raw_ci,)
        if any(ci is None for ci in ci_vectors):
            raise ValueError("PySCF CASSCF object does not contain CI roots.")
        raw_energies = getattr(mc, "e_states", None)
        if raw_energies is None:
            raw_energies = (getattr(mc, "e_tot"),)
        energies = np.atleast_1d(np.asarray(raw_energies, dtype=float))
        if len(energies) != len(ci_vectors) or not np.all(np.isfinite(energies)):
            raise ValueError("PySCF root energies and CI-vector counts differ.")
        nroot = len(ci_vectors)
        labels = (
            tuple(f"R{index + 1}" for index in range(nroot))
            if root_labels is None
            else tuple(str(label) for label in root_labels)
        )
        if len(labels) != nroot or any(not label.strip() for label in labels):
            raise ValueError("one nonempty label is required per PySCF root.")
        ncore = int(mc.ncore)
        ncas = int(mc.ncas)
        nelecas = tuple(int(value) for value in mc.nelecas)
        nmo = int(np.asarray(mc.mo_coeff).shape[1])
        reference_ms_twice = nelecas[0] - nelecas[1]
        spin_data = _root_spin_data_v241(
            ci_vectors,
            ncas=ncas,
            nelecas=nelecas,
            declared=root_spin_twice,
        )
        roots = tuple(
            SpinFreeRootV241(
                label=label,
                energy_hartree=float(energy),
                spin_twice=spin_twice,
                reference_ms_twice=reference_ms_twice,
                spin_square=spin_square,
            ).validate()
            for label, energy, (spin_twice, spin_square) in zip(
                labels, energies, spin_data
            )
        )
        microstates = complete_spin_microstates_v241(roots)
        state_weights = _normalised_weights_v241(
            getattr(mc, "weights", None) if weights is None else weights,
            nroot,
        )
        density_mo = state_average_density_mo_from_pyscf_ci_v241(
            ci_vectors,
            ncore=ncore,
            ncas=ncas,
            nmo=nmo,
            nelecas=nelecas,
            weights=state_weights,
        )
        integrals = build_pyscf_bp_somf_integrals_v241(
            self.mol, mc.mo_coeff, density_mo
        )
        reduced, sample_ms, reference_coefficients = (
            wigner_reduced_transition_density_from_pyscf_ci_v241(
                ci_vectors,
                roots,
                ncore=ncore,
                ncas=ncas,
                nmo=nmo,
                nelecas=nelecas,
            )
        )
        matrices = assemble_state_interaction_soc_v241(
            roots, reduced, integrals.effective_mo_cartesian
        )

        atomic_symbols = tuple(self.mol.atom_symbol(i) for i in range(self.mol.natm))
        if isotope_masses_amu is None:
            isotope_masses_amu = tuple(
                float(value) for value in self.mol.atom_mass_list()
            )
        else:
            isotope_masses_amu = tuple(float(value) for value in isotope_masses_amu)
        basis_label = str(self.mol.basis) if basis_label is None else str(basis_label)
        calculation_input = {
            "atom_symbols": list(atomic_symbols),
            "geometry_bohr": np.asarray(self.mol.atom_coords(), dtype=float).tolist(),
            "charge": int(self.mol.charge),
            "spin_twice": int(self.mol.spin),
            "basis": basis_label,
            "ncore": ncore,
            "ncas": ncas,
            "nelecas": list(nelecas),
            "weights": state_weights.tolist(),
            "root_labels": list(labels),
            "root_spin_twice": [root.spin_twice for root in roots],
            "operator": BP_SOMF_OPERATOR_FAMILY_V241,
            "scalar_relativistic_method": "none",
        }
        calculation_input_sha256 = _sha256_v241(calculation_input)
        capabilities = MolecularSOCCapabilitiesV230(
            static_soc=True,
            spin_free_derivatives=False,
            soc_derivatives=False,
            derivative_connections=False,
            cross_geometry_overlaps=False,
            deterministic_replay=False,
            analytic_soc_derivatives=False,
        ).validate()
        identity = MolecularSOCBackendIdentityV230(
            backend_name="PySCF",
            backend_version=PYSCF_REQUIRED_VERSION_V241,
            source_kind="live_ab_initio",
            electronic_method=(
                f"common-orbital SA-CASSCF({sum(nelecas)}e,{ncas}o) state interaction"
            ),
            basis=basis_label,
            charge=int(self.mol.charge),
            electron_count=int(self.mol.nelectron),
            soc_operator=BP_SOMF_OPERATOR_FAMILY_V241,
            scalar_relativistic_method="none (nonrelativistic scalar Hamiltonian)",
            derivative_method=BP_SOMF_STATIC_LIMITATION_V241,
            active_space=f"CAS({sum(nelecas)}e,{ncas}o)",
            molecule_name=str(molecule_name),
            atom_symbols=atomic_symbols,
            isotope_masses_amu=isotope_masses_amu,
            reference_geometry_bohr=tuple(
                tuple(float(value) for value in row)
                for row in np.asarray(self.mol.atom_coords(), dtype=float)
            ),
            calculation_input_sha256=calculation_input_sha256,
            environment_sha256=environment_sha256,
            extra={
                "provider": self.adapter_name,
                "provider_version": self.adapter_version,
                "state_average_weights": state_weights.tolist(),
                "root_spin_twice": [root.spin_twice for root in roots],
                "complete_microstate_order": list(matrices.state_order),
                "one_electron_integral": BP_SOMF_ONE_ELECTRON_INTEGRAL_V241,
                "two_electron_integral": BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
                "static_limitation": BP_SOMF_STATIC_LIMITATION_V241,
            },
        ).validate()
        molecular_contract = MolecularSOCAdmissionContractV230(
            capabilities=capabilities,
            identity=identity,
            state_tracking_policy="not available for static-only v0.24.1 result",
            coordinate_definition="Cartesian molecular geometry in bohr",
            all_electronic_calculations_converged=True,
        ).validate()
        parity = identity.electron_parity
        convention = MolecularSOCMatrixConventionV233(
            schema=SOC_CONVENTION_SCHEMA_V233,
            operator_family=BP_SOMF_OPERATOR_FAMILY_V241,
            one_electron_treatment=(
                "all-electron nuclear-attraction SOC from PySCF int1e_prinvxp"
            ),
            two_electron_treatment=(
                "PySCF int2e_p1vxp1 contracted as J-3/2 K_left-3/2 K_right"
            ),
            mean_field_approximation="state-averaged spin-free one-particle density",
            prefactor_convention=(
                "0.5/c^2 included once; real antisymmetric Cartesian operator "
                "multiplied by -i in the scalar-MO basis"
            ),
            scalar_relativistic_method=identity.scalar_relativistic_method,
            source_basis="common real scalar SA-CASSCF molecular orbitals",
            target_basis="complete |root,S,M_S> spin-microstate basis",
            state_order=matrices.state_order,
            electron_parity=parity,
        ).validate()
        symmetry = SOCSymmetryContractV221(
            parity,
            matrices.time_reversal_matrix,
            matrices.root_projectors,
            external_magnetic_field=False,
        )
        model_space = _model_space_v241(microstates, self.mol.charge)
        provenance = ElectronicOperatorProvenanceV213(
            model_name="v0.24.1 PySCF direct static molecular SOC",
            model_version="1",
            model_space=model_space,
            spin_free_method=identity.electronic_method,
            soc_enabled=True,
            soc_method=BP_SOMF_OPERATOR_FAMILY_V241,
            scalar_relativistic_method=identity.scalar_relativistic_method,
            derivative_method=BP_SOMF_STATIC_LIMITATION_V241,
            parameters={
                **symmetry.as_provenance_parameters(),
                "backend_identity": identity.as_dict(),
                "capabilities": capabilities.as_dict(),
                "soc_convention": convention.as_dict(),
                "soc_convention_fingerprint": convention.fingerprint(),
                "calculation_input": calculation_input,
                "calculation_input_sha256": calculation_input_sha256,
                "nac_convention": (
                    "not exercised by static SOC; inherited dynamics convention unchanged"
                ),
            },
        ).validate()
        self.capabilities = capabilities
        self.identity = identity
        self._molecular_soc_contract = molecular_contract
        self._soc_symmetry_contract = symmetry
        self.convention = convention
        self.provenance = provenance
        self.result = PySCFStateInteractionSOCResultV241(
            matrices=matrices,
            integrals=integrals,
            wigner_reduced_transition_density=reduced,
            wigner_sample_ms_twice=sample_ms,
            wigner_reference_coefficients=reference_coefficients,
            capabilities=capabilities,
            identity=identity,
            molecular_soc_contract=molecular_contract,
            convention=convention,
            provenance=provenance,
            scf_converged=True,
            casscf_converged=True,
            soc_assembled=True,
            runtime_probe=self.runtime_probe,
            metadata={
                "calculation_input_sha256": calculation_input_sha256,
                "state_average_weights": state_weights.tolist(),
                "ncore": ncore,
                "ncas": ncas,
                "nmo": nmo,
                "nelecas": list(nelecas),
            },
        ).validate()

    @property
    def backend_version(self):
        return self.runtime_probe.module_version

    @property
    def molecular_soc_contract(self):
        return self._molecular_soc_contract

    @property
    def soc_symmetry_contract(self):
        return self._soc_symmetry_contract

    @property
    def time_reversal_matrix(self):
        return self.result.matrices.time_reversal_matrix.copy()

    @property
    def projectors(self):
        return {
            name: matrix.copy()
            for name, matrix in self.result.matrices.root_projectors.items()
        }

    def evaluate_static_soc(self):
        return self.result.validate()

    def components(self, q):
        raise RuntimeError(
            "v0.24.1 PySCF BP-SOMF is static-only; physical H_SOC derivatives "
            "are required before components(q) can be exposed."
        )

    def evaluate_snapshot(self, q):
        raise RuntimeError(
            "v0.24.1 PySCF BP-SOMF is not trajectory-ready; use "
            "evaluate_static_soc() for the fixed-geometry matrix."
        )

    def snapshot_overlap(self, left, right):
        raise RuntimeError(
            "v0.24.1 PySCF BP-SOMF does not implement cross-geometry "
            "many-electron overlaps."
        )


def tampered_soc_result_v241(result, **changes):
    """Test helper that preserves frozen dataclass semantics while changing fields."""

    if not isinstance(result, PySCFStateInteractionSOCResultV241):
        raise TypeError("tampering helper requires a v0.24.1 PySCF SOC result.")
    return replace(result, **changes)
