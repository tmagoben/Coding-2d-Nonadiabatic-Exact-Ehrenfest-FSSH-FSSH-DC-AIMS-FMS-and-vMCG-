"""Frozen OpenMolcas RASSI-SO snapshot protocol for v0.24.0.

The protocol is intentionally narrower than a general OpenMolcas interface.  It
defines the first molecular-SOC intake target and is suitable as an out-of-band
trust anchor.  It does not claim that an OpenMolcas calculation has been run.
"""

from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np

from .analytic_soc_models_v220 import (
    singlet_triplet_projectors_v220,
    singlet_triplet_time_reversal_matrix_v220,
)
from .electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicStateDescriptorV213,
)
from .molecular_soc_convention_v233 import (
    MolecularSOCMatrixConventionV233,
    SOC_CONVENTION_SCHEMA_V233,
)
from .soc_admission_v221 import SOCSymmetryContractV221


OPENMOLCAS_PROTOCOL_SCHEMA_V240 = "gnd-openmolcas-rassi-so-protocol-v0.24.0"
OPENMOLCAS_EXPORT_SCHEMA_V240 = "gnd-openmolcas-rassi-so-export-v0.24.0"
OPENMOLCAS_MANIFEST_SCHEMA_V240 = "gnd-openmolcas-rassi-so-bundle-v0.24.0"
OPENMOLCAS_ADAPTER_NAME_V240 = "gnd-openmolcas-rassi-snapshot"
OPENMOLCAS_ADAPTER_VERSION_V240 = "0.24.0"


def _canonical_bytes_v240(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class OpenMolcasRASSIProtocolV240:
    """Exact method and artifact contract for the first external SOC intake."""

    schema: str
    backend_name: str
    backend_version: str
    molecule_name: str
    atom_symbols: tuple[str, ...]
    isotope_masses_amu: tuple[float, ...]
    reference_geometry_bohr: tuple[tuple[float, float, float], ...]
    charge: int
    electron_count: int
    point_group: str
    basis: str
    scalar_relativistic_method: str
    active_space: str
    spin_free_method: str
    dynamic_correlation_method: str
    soc_operator: str
    rassi_state_selection: str
    rassi_energy_source: str
    state_order: tuple[str, ...]
    displacement_steps_bohr: tuple[float, ...]
    cross_geometry_overlap_method: str
    geometry_unit: str = "bohr"
    energy_unit: str = "hartree"
    derivative_unit: str = "hartree/bohr"
    spin_quantization_axis: str = "z"
    external_magnetic_field: bool = False

    @property
    def coordinate_dimension(self):
        return 3 * len(self.atom_symbols)

    def validate(self):
        if self.schema != OPENMOLCAS_PROTOCOL_SCHEMA_V240:
            raise ValueError("OpenMolcas protocol schema mismatch.")
        for name in (
            "backend_name",
            "backend_version",
            "molecule_name",
            "point_group",
            "basis",
            "scalar_relativistic_method",
            "active_space",
            "spin_free_method",
            "dynamic_correlation_method",
            "soc_operator",
            "rassi_state_selection",
            "rassi_energy_source",
            "cross_geometry_overlap_method",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"OpenMolcas protocol {name} cannot be empty.")
        if self.backend_name != "OpenMolcas" or self.backend_version != "26.06":
            raise ValueError("v0.24.0 pins OpenMolcas 26.06 exactly.")
        if int(self.charge) != self.charge:
            raise ValueError("molecular charge must be an integer.")
        if int(self.electron_count) != self.electron_count or self.electron_count < 1:
            raise ValueError("electron_count must be a positive integer.")
        if self.electron_count % 2:
            raise ValueError("the first v0.24.0 target is the even-electron sector.")
        symbols = tuple(str(item).strip() for item in self.atom_symbols)
        masses = np.asarray(self.isotope_masses_amu, dtype=float)
        geometry = np.asarray(self.reference_geometry_bohr, dtype=float)
        if not symbols or any(not item for item in symbols):
            raise ValueError("atom_symbols must be nonempty.")
        if masses.shape != (len(symbols),) or np.any(masses <= 0.0):
            raise ValueError("one positive isotope mass is required per atom.")
        if geometry.shape != (len(symbols), 3) or not np.all(np.isfinite(geometry)):
            raise ValueError("reference geometry must be finite with shape (natom,3).")
        if self.geometry_unit != "bohr" or self.energy_unit != "hartree":
            raise ValueError("the frozen internal units are bohr and hartree.")
        if self.derivative_unit != "hartree/bohr":
            raise ValueError("the frozen derivative unit is hartree/bohr.")
        if self.spin_quantization_axis != "z":
            raise ValueError("v0.24.0 freezes the spin quantization axis to z.")
        if type(self.external_magnetic_field) is not bool:
            raise TypeError("external_magnetic_field must be a native Boolean.")
        if self.external_magnetic_field:
            raise ValueError("the first external snapshot excludes magnetic fields.")
        if tuple(self.state_order) != (
            "S0(M=0)",
            "T1(M=-1)",
            "T1(M=0)",
            "T1(M=+1)",
        ):
            raise ValueError("v0.24.0 requires one singlet and one complete triplet.")
        steps = np.asarray(self.displacement_steps_bohr, dtype=float)
        if steps.shape != (3,) or not np.all(np.isfinite(steps)):
            raise ValueError("exactly three finite displacement steps are required.")
        if np.any(steps <= 0.0) or not np.all(np.diff(steps) < 0.0):
            raise ValueError("displacement steps must be positive and strictly decreasing.")
        return self

    def as_dict(self):
        self.validate()
        payload = asdict(self)
        for name in (
            "atom_symbols",
            "isotope_masses_amu",
            "reference_geometry_bohr",
            "state_order",
            "displacement_steps_bohr",
        ):
            payload[name] = np.asarray(payload[name]).tolist()
        return payload

    def fingerprint(self):
        return hashlib.sha256(_canonical_bytes_v240(self.as_dict())).hexdigest()

    def expected_record_ids(self):
        identifiers = ["reference"]
        for coordinate in range(self.coordinate_dimension):
            for step_index in range(len(self.displacement_steps_bohr)):
                identifiers.extend(
                    (
                        f"q{coordinate:02d}_h{step_index}_minus",
                        f"q{coordinate:02d}_h{step_index}_plus",
                    )
                )
        return tuple(identifiers)

    def model_space(self):
        return ElectronicModelSpaceV213(
            name="v0.24.0 OpenMolcas water S0/T1 intake space",
            representation="fixed_spin_diabatic",
            states=(
                ElectronicStateDescriptorV213("S0(M=0)", "S0", 1, "M=0", 0),
                ElectronicStateDescriptorV213("T1(M=-1)", "T1", 3, "M=-1", 0),
                ElectronicStateDescriptorV213("T1(M=0)", "T1", 3, "M=0", 0),
                ElectronicStateDescriptorV213("T1(M=+1)", "T1", 3, "M=+1", 0),
            ),
            complete_multiplets=True,
        ).validate()

    def symmetry_contract(self):
        return SOCSymmetryContractV221(
            electron_parity="even",
            time_reversal_matrix=singlet_triplet_time_reversal_matrix_v220(),
            projectors=singlet_triplet_projectors_v220(),
        )

    def soc_convention(self):
        return MolecularSOCMatrixConventionV233(
            schema=SOC_CONVENTION_SCHEMA_V233,
            operator_family=self.soc_operator,
            one_electron_treatment="OpenMolcas SEWARD AMFI integrals",
            two_electron_treatment="AMFI atomic mean-field approximation",
            mean_field_approximation="atomic mean-field integrals (AMFI)",
            prefactor_convention=(
                "OpenMolcas RASSI SPINORBIT matrix elements in hartree; "
                "no adapter-side prefactor"
            ),
            scalar_relativistic_method=self.scalar_relativistic_method,
            source_basis="RASSI orthonormal spin-free eigenstates",
            target_basis="fixed ordered S0/T1 spin-component frame",
            state_order=self.state_order,
            electron_parity="even",
        ).validate()


def water_rassi_so_protocol_v240():
    """Return the sole v0.24.0 method trust anchor.

    The geometry is a reproducible near-equilibrium water geometry.  Method choices
    are project choices; their presence here is not evidence of execution or accuracy.
    """

    return OpenMolcasRASSIProtocolV240(
        schema=OPENMOLCAS_PROTOCOL_SCHEMA_V240,
        backend_name="OpenMolcas",
        backend_version="26.06",
        molecule_name="water",
        atom_symbols=("O", "H", "H"),
        isotope_masses_amu=(15.99491461957, 1.00782503223, 1.00782503223),
        reference_geometry_bohr=(
            (0.0, 0.0, 0.0),
            (0.0, 1.430824928, 1.107688329),
            (0.0, -1.430824928, 1.107688329),
        ),
        charge=0,
        electron_count=10,
        point_group="C1 for every reference and displaced calculation",
        basis="ANO-RCC-VDZP",
        scalar_relativistic_method="second-order Douglas-Kroll-Hess (DKH2)",
        active_space="CAS(8,6), identical orbital/state definition at every geometry",
        spin_free_method="state-specific CASSCF followed by single-state CASPT2",
        dynamic_correlation_method="single-state CASPT2 diagonal energies via EJOB",
        soc_operator="OpenMolcas RASSI-SO with SEWARD AMFI",
        rassi_state_selection="NROFJOBIPH=2 1 1;1;1",
        rassi_energy_source="EJOB",
        state_order=("S0(M=0)", "T1(M=-1)", "T1(M=0)", "T1(M=+1)"),
        displacement_steps_bohr=(0.004, 0.002, 0.001),
        cross_geometry_overlap_method=(
            "independent biorthogonal CASSCF wavefunction overlap exporter; "
            "not inferred from the within-geometry RASSI overlap matrix"
        ),
    ).validate()


def openmolcas_protocol_from_dict_v240(payload):
    if not isinstance(payload, dict):
        raise TypeError("OpenMolcas protocol payload must be a mapping.")
    expected = set(OpenMolcasRASSIProtocolV240.__dataclass_fields__)
    if set(payload) != expected:
        raise ValueError("OpenMolcas protocol field set mismatch.")
    converted = dict(payload)
    for name in (
        "atom_symbols",
        "isotope_masses_amu",
        "reference_geometry_bohr",
        "state_order",
        "displacement_steps_bohr",
    ):
        converted[name] = tuple(converted[name])
    converted["reference_geometry_bohr"] = tuple(
        tuple(row) for row in converted["reference_geometry_bohr"]
    )
    return OpenMolcasRASSIProtocolV240(**converted).validate()
