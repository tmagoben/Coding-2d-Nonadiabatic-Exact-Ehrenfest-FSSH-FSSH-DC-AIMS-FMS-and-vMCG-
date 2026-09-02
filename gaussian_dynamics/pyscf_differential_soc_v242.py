"""Connected-geometry PySCF BP-SOMF differentials for v0.24.2.

This module advances the v0.24.1 fixed-geometry implementation in two deliberately
separate ways:

* the production SOMF contraction uses PySCF's direct JK driver and never materializes
  the rank-five two-electron SOC tensor;
* neighboring complete spin-microstate spaces are connected by exact restricted-
  CASSCF many-electron overlaps and their certified unitary polar transports before
  any finite difference is formed.

The result is differential evidence, not yet a trajectory-ready provider.  Analytic
SOC derivatives, a full Cartesian molecular scan, and accuracy admission remain
outside the v0.24.2 capability boundary.
"""

from dataclasses import dataclass, field
import hashlib
import json
import math

import numpy as np

from .finite_manifold_transport_v233 import (
    FiniteManifoldOverlapPolicyV233,
    certified_transport_from_overlap_v233,
)
from .pyscf_runtime_v232 import guarded_pyscf_runtime_v232
from .pyscf_state_interaction_soc_v241 import (
    BP_SOMF_CARTESIAN_ORDER_V241,
    BP_SOMF_ONE_ELECTRON_INTEGRAL_V241,
    BP_SOMF_OPERATOR_FAMILY_V241,
    BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    BPSOMFIntegralsV241,
    SpinFreeRootV241,
    assemble_state_interaction_soc_v241,
    build_pyscf_bp_somf_integrals_v241,
    complete_spin_microstates_v241,
    require_pyscf_static_soc_runtime_v241,
    state_average_density_mo_from_pyscf_ci_v241,
    wigner_reduced_transition_density_from_pyscf_ci_v241,
)
from .pyscf_wavefunction_overlap import (
    CASSCFWavefunctionSnapshot,
    casscf_state_overlap_matrix,
)


PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242 = (
    "gnd-pyscf-connected-geometry-soc-differential-v0.24.2"
)
PYSCF_DIRECT_JK_SOMF_STRATEGY_V242 = (
    "PySCF direct JK int2e_p1vxp1 contraction; no rank-five tensor materialization"
)
PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242 = "connected_geometry_differential_preview"
OH_BOND_LENGTH_BOHR_V242 = 1.83256418024373
OH_ISOTOPE_MASSES_AMU_V242 = (15.99491461957, 1.00782503223)
# The 0.08 -> 0.04 -> 0.02 bohr ladder lies on the observed centered-difference
# truncation plateau for both H_sf and H_SOC.  Finer steps enter independent CASSCF
# convergence noise and are deliberately not presented as improved evidence.
OH_BOND_STEPS_BOHR_V242 = (8.0e-2, 4.0e-2, 2.0e-2)


def _canonical_v242(value):
    if isinstance(value, np.generic):
        return _canonical_v242(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_v242(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_v242(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v242(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("canonical v0.24.2 evidence cannot contain non-finite data.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.24.2 canonical value: {type(value).__name__}")


def _canonical_bytes_v242(value):
    return json.dumps(
        _canonical_v242(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def _sha256_v242(value):
    return hashlib.sha256(_canonical_bytes_v242(value)).hexdigest()


def _array_sha256_v242(*arrays):
    digest = hashlib.sha256()
    for array in arrays:
        item = np.ascontiguousarray(np.asarray(array))
        digest.update(item.dtype.str.encode("ascii"))
        digest.update(str(tuple(item.shape)).encode("ascii"))
        digest.update(item.tobytes(order="C"))
    return digest.hexdigest()


def _complex_pairs_v242(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_frobenius_v242(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right, ord="fro") / scale)


def _normalised_weights_v242(weights, nroot):
    if weights is None:
        result = np.full(int(nroot), 1.0 / int(nroot), dtype=float)
    else:
        result = np.asarray(weights, dtype=float)
    if result.shape != (int(nroot),):
        raise ValueError("one state-average weight is required per root.")
    if not np.all(np.isfinite(result)) or np.any(result < 0.0):
        raise ValueError("state-average weights must be finite and nonnegative.")
    total = float(np.sum(result))
    if total <= 0.0:
        raise ValueError("state-average weights must have positive total weight.")
    result = result / total
    if abs(float(np.sum(result)) - 1.0) > 1.0e-14:
        raise ValueError("normalized state-average weights do not sum to one.")
    return result


def build_pyscf_bp_somf_integrals_direct_jk_v242(
    mol,
    mo_coeff,
    state_average_density_mo,
    *,
    imaginary_tolerance=1.0e-12,
):
    """Build BP-SOMF integrals through the direct PySCF JK contraction path."""

    require_pyscf_static_soc_runtime_v241()
    from pyscf.lib.parameters import LIGHT_SPEED
    from pyscf.scf import jk

    coefficients = np.asarray(mo_coeff, dtype=complex)
    density_mo = np.asarray(state_average_density_mo, dtype=complex)
    nao = int(mol.nao_nr())
    if coefficients.ndim != 2 or coefficients.shape[0] != nao:
        raise ValueError("MO coefficients have the wrong AO dimension.")
    if density_mo.shape != (coefficients.shape[1], coefficients.shape[1]):
        raise ValueError("state-average MO density has the wrong shape.")
    if not np.all(np.isfinite(coefficients)) or not np.all(np.isfinite(density_mo)):
        raise ValueError("MO coefficients and density must be finite.")
    if np.max(np.abs(coefficients.imag)) > float(imaginary_tolerance):
        raise ValueError("v0.24.2 direct-JK SOMF requires real common scalar orbitals.")
    if np.max(np.abs(density_mo.imag)) > float(imaginary_tolerance):
        raise ValueError("v0.24.2 direct-JK SOMF requires a real spin-free density.")
    coefficients = np.asarray(coefficients.real, dtype=float)
    density_mo = np.asarray(density_mo.real, dtype=float)
    density_ao = coefficients @ density_mo @ coefficients.T

    one_electron = np.zeros((3, nao, nao), dtype=float)
    for atom_index in range(int(mol.natm)):
        with mol.with_rinv_origin(mol.atom_coord(atom_index)):
            one_electron += float(mol.atom_charge(atom_index)) * np.asarray(
                mol.intor(BP_SOMF_ONE_ELECTRON_INTEGRAL_V241, comp=3),
                dtype=float,
            )

    coulomb_like, exchange_left, exchange_right = jk.get_jk(
        mol,
        [density_ao, density_ao, density_ao],
        scripts=(
            "ijkl,kl->ij",
            "ijkl,jk->il",
            "ijkl,li->kj",
        ),
        intor=BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    )
    two_electron_somf = (
        np.asarray(coulomb_like, dtype=float)
        - 1.5 * np.asarray(exchange_left, dtype=float)
        - 1.5 * np.asarray(exchange_right, dtype=float)
    )
    if two_electron_somf.shape != (3, nao, nao):
        raise ValueError("PySCF direct JK returned an unexpected SOC tensor shape.")

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
        one_electron_ao_cartesian=one_electron,
        two_electron_somf_ao_cartesian=two_electron_somf,
        effective_ao_cartesian=effective_ao,
        effective_mo_cartesian=effective_mo,
        state_average_density_mo=density_mo,
        state_average_density_ao=density_ao,
        light_speed_au=light_speed,
        prefactor=prefactor,
        cartesian_component_order=BP_SOMF_CARTESIAN_ORDER_V241,
        one_electron_integral=BP_SOMF_ONE_ELECTRON_INTEGRAL_V241,
        two_electron_integral=BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    ).validate()


def _spin_roots_from_casscf_v242(
    ci_vectors,
    energies,
    *,
    ncas,
    nelecas,
    root_labels,
    root_spin_twice=None,
):
    from pyscf.fci import spin_op

    labels = tuple(str(label) for label in root_labels)
    if len(labels) != len(ci_vectors) or len(set(labels)) != len(labels):
        raise ValueError("unique root labels are required for every CI root.")
    declared = None if root_spin_twice is None else tuple(int(x) for x in root_spin_twice)
    if declared is not None and len(declared) != len(ci_vectors):
        raise ValueError("one declared spin is required per CI root.")
    reference_ms_twice = int(nelecas[0]) - int(nelecas[1])
    roots = []
    for index, (label, energy, ci) in enumerate(zip(labels, energies, ci_vectors)):
        spin_square, multiplicity = spin_op.spin_square(ci, int(ncas), nelecas)
        observed_spin_twice = int(round(float(multiplicity) - 1.0))
        spin_twice = observed_spin_twice if declared is None else declared[index]
        if declared is not None and spin_twice != observed_spin_twice:
            raise ValueError("declared root spin disagrees with PySCF <S^2>.")
        roots.append(
            SpinFreeRootV241(
                label=label,
                energy_hartree=float(energy),
                spin_twice=spin_twice,
                reference_ms_twice=reference_ms_twice,
                spin_square=float(spin_square),
            ).validate()
        )
    return tuple(roots)


@dataclass(frozen=True)
class PySCFSOCGeometrySnapshotV242:
    """One converged geometry with direct-JK SOC and overlap-capable CI data."""

    geometry_bohr: np.ndarray
    roots: tuple
    matrices: object
    integrals: BPSOMFIntegralsV241
    wigner_reduced_transition_density: np.ndarray
    wigner_sample_ms_twice: np.ndarray
    wigner_reference_coefficients: np.ndarray
    wavefunction_snapshot: CASSCFWavefunctionSnapshot
    state_average_weights: np.ndarray
    environment_sha256: str
    calculation_input: dict
    scf_energy_hartree: float
    casscf_object: object = field(repr=False, compare=False)

    def __post_init__(self):
        object.__setattr__(self, "geometry_bohr", np.asarray(self.geometry_bohr, dtype=float).copy())
        object.__setattr__(
            self,
            "wigner_reduced_transition_density",
            np.asarray(self.wigner_reduced_transition_density, dtype=complex).copy(),
        )
        object.__setattr__(
            self,
            "wigner_sample_ms_twice",
            np.asarray(self.wigner_sample_ms_twice, dtype=int).copy(),
        )
        object.__setattr__(
            self,
            "wigner_reference_coefficients",
            np.asarray(self.wigner_reference_coefficients, dtype=float).copy(),
        )
        object.__setattr__(
            self, "state_average_weights", np.asarray(self.state_average_weights, dtype=float).copy()
        )
        object.__setattr__(self, "roots", tuple(self.roots))
        object.__setattr__(self, "calculation_input", dict(self.calculation_input))

    @property
    def state_order(self):
        return self.matrices.state_order

    @property
    def wavefunction_sha256(self):
        arrays = [self.wavefunction_snapshot.mo_coeff]
        arrays.extend(self.wavefunction_snapshot.ci_roots)
        return _array_sha256_v242(*arrays)

    def validate(self, tolerance=1.0e-10):
        geometry = np.asarray(self.geometry_bohr, dtype=float)
        molecular_geometry = np.asarray(self.wavefunction_snapshot.mol.atom_coords(), dtype=float)
        if geometry.ndim != 2 or geometry.shape[1] != 3 or not np.all(np.isfinite(geometry)):
            raise ValueError("snapshot geometry must be a finite (natom,3) array.")
        if geometry.shape != molecular_geometry.shape or np.max(np.abs(geometry - molecular_geometry)) > 1.0e-12:
            raise ValueError("snapshot geometry disagrees with its PySCF molecule.")
        roots = tuple(root.validate() for root in self.roots)
        if not roots or tuple(self.matrices.roots) != roots:
            raise ValueError("snapshot roots disagree with state-interaction matrices.")
        self.matrices.validate(tolerance=tolerance)
        self.integrals.validate(tolerance=tolerance)
        if self.wavefunction_snapshot.nroots != len(roots):
            raise ValueError("wavefunction and SOC root counts differ.")
        if tuple(self.state_order) != tuple(state.label for state in complete_spin_microstates_v241(roots)):
            raise ValueError("snapshot complete-multiplet order is inconsistent.")
        weights = np.asarray(self.state_average_weights, dtype=float)
        if weights.shape != (len(roots),) or abs(float(np.sum(weights)) - 1.0) > 1.0e-12:
            raise ValueError("snapshot state-average weights are invalid.")
        if not isinstance(self.environment_sha256, str) or len(self.environment_sha256) != 64:
            raise ValueError("snapshot environment identity must be a SHA-256 digest.")
        if any(character not in "0123456789abcdef" for character in self.environment_sha256):
            raise ValueError("snapshot environment identity must be lowercase hexadecimal.")
        if self.calculation_input.get("somf_contraction") != PYSCF_DIRECT_JK_SOMF_STRATEGY_V242:
            raise ValueError("snapshot does not declare the v0.24.2 direct-JK strategy.")
        if not np.isfinite(float(self.scf_energy_hartree)):
            raise ValueError("snapshot SCF energy must be finite.")
        return self

    def compact_dict(self):
        self.validate()
        return {
            "schema": PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242,
            "geometry_bohr": self.geometry_bohr.tolist(),
            "roots": [root.as_dict() for root in self.roots],
            "state_order": list(self.state_order),
            "state_average_weights": self.state_average_weights.tolist(),
            "environment_sha256": self.environment_sha256,
            "calculation_input": self.calculation_input,
            "scf_energy_hartree": float(self.scf_energy_hartree),
            "wavefunction_sha256": self.wavefunction_sha256,
            "integrals": self.integrals.as_dict(include_matrices=False),
            "matrices": self.matrices.as_dict(include_matrices=True),
            "wigner_sample_ms_twice": self.wigner_sample_ms_twice.tolist(),
            "wigner_reference_coefficients": self.wigner_reference_coefficients.tolist(),
        }

    def fingerprint(self):
        return _sha256_v242(self.compact_dict())


def build_pyscf_soc_geometry_snapshot_v242(
    casscf,
    *,
    environment_sha256,
    root_labels=None,
    root_spin_twice=None,
    weights=None,
    molecule_name="PySCF molecule",
    basis_label=None,
):
    """Build one overlap-capable direct-JK state-interaction SOC snapshot."""

    require_pyscf_static_soc_runtime_v241()
    if not bool(getattr(getattr(casscf, "_scf", None), "converged", False)):
        raise RuntimeError("PySCF SCF reference is not converged.")
    if not bool(getattr(casscf, "converged", False)):
        raise RuntimeError("PySCF CASSCF calculation is not converged.")
    mol = getattr(casscf, "mol", None)
    if mol is None:
        raise TypeError("CASSCF snapshot requires a PySCF molecule.")
    raw_ci = getattr(casscf, "ci", None)
    ci_vectors = tuple(raw_ci) if isinstance(raw_ci, (tuple, list)) else (raw_ci,)
    if not ci_vectors or any(ci is None for ci in ci_vectors):
        raise ValueError("CASSCF snapshot does not contain CI roots.")
    ci_vectors = tuple(np.asarray(ci) for ci in ci_vectors)
    raw_energies = getattr(casscf, "e_states", None)
    if raw_energies is None:
        raw_energies = (getattr(casscf, "e_tot"),)
    energies = np.atleast_1d(np.asarray(raw_energies, dtype=float))
    if energies.shape != (len(ci_vectors),) or not np.all(np.isfinite(energies)):
        raise ValueError("CASSCF state energies and CI roots are inconsistent.")
    nroot = len(ci_vectors)
    labels = (
        tuple(f"R{index + 1}" for index in range(nroot))
        if root_labels is None
        else tuple(str(label) for label in root_labels)
    )
    ncore = int(casscf.ncore)
    ncas = int(casscf.ncas)
    nelecas = tuple(int(value) for value in casscf.nelecas)
    nmo = int(np.asarray(casscf.mo_coeff).shape[1])
    roots = _spin_roots_from_casscf_v242(
        ci_vectors,
        energies,
        ncas=ncas,
        nelecas=nelecas,
        root_labels=labels,
        root_spin_twice=root_spin_twice,
    )
    state_weights = _normalised_weights_v242(
        getattr(casscf, "weights", None) if weights is None else weights,
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
    integrals = build_pyscf_bp_somf_integrals_direct_jk_v242(
        mol, casscf.mo_coeff, density_mo
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
    geometry = np.asarray(mol.atom_coords(), dtype=float)
    calculation_input = {
        "molecule_name": str(molecule_name),
        "atom_symbols": [mol.atom_symbol(index) for index in range(int(mol.natm))],
        "geometry_bohr": geometry.tolist(),
        "charge": int(mol.charge),
        "spin_twice": int(mol.spin),
        "basis": str(mol.basis) if basis_label is None else str(basis_label),
        "ncore": ncore,
        "ncas": ncas,
        "nelecas": list(nelecas),
        "state_average_weights": state_weights.tolist(),
        "root_labels": list(labels),
        "root_spin_twice": [root.spin_twice for root in roots],
        "operator": BP_SOMF_OPERATOR_FAMILY_V241,
        "somf_contraction": PYSCF_DIRECT_JK_SOMF_STRATEGY_V242,
        "scalar_relativistic_method": "none",
    }
    wavefunction = CASSCFWavefunctionSnapshot(
        mol=mol,
        mo_coeff=np.asarray(casscf.mo_coeff).copy(),
        ci_roots=tuple(np.asarray(ci).copy() for ci in ci_vectors),
        ncore=ncore,
        ncas=ncas,
        nelecas=nelecas,
        metadata={
            "schema": PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242,
            "calculation_input_sha256": _sha256_v242(calculation_input),
            "root_labels": list(labels),
        },
    )
    return PySCFSOCGeometrySnapshotV242(
        geometry_bohr=geometry,
        roots=roots,
        matrices=matrices,
        integrals=integrals,
        wigner_reduced_transition_density=reduced,
        wigner_sample_ms_twice=sample_ms,
        wigner_reference_coefficients=reference_coefficients,
        wavefunction_snapshot=wavefunction,
        state_average_weights=state_weights,
        environment_sha256=str(environment_sha256),
        calculation_input=calculation_input,
        scf_energy_hartree=float(casscf._scf.e_tot),
        casscf_object=casscf,
    ).validate()


def complete_multiplet_overlap_v242(left, right):
    """Lift a spin-free CASSCF root overlap into complete |root,S,M_S> spaces."""

    left = left.validate()
    right = right.validate()
    if len(left.roots) != len(right.roots):
        raise ValueError("endpoint root counts differ.")
    if tuple(root.label for root in left.roots) != tuple(root.label for root in right.roots):
        raise ValueError("endpoint root labels differ.")
    if tuple(root.spin_twice for root in left.roots) != tuple(
        root.spin_twice for root in right.roots
    ):
        raise ValueError("endpoint root-spin declarations differ.")
    root_overlap = np.asarray(
        casscf_state_overlap_matrix(
            left.wavefunction_snapshot, right.wavefunction_snapshot
        ),
        dtype=complex,
    )
    left_states = left.matrices.microstates
    right_states = right.matrices.microstates
    result = np.zeros((len(left_states), len(right_states)), dtype=complex)
    for row, bra in enumerate(left_states):
        for column, ket in enumerate(right_states):
            if bra.spin_twice == ket.spin_twice and bra.ms_twice == ket.ms_twice:
                result[row, column] = root_overlap[bra.root_index, ket.root_index]
    if not np.all(np.isfinite(result)):
        raise ValueError("complete-multiplet overlap contains non-finite data.")
    return result


def phase_align_complete_multiplet_overlap_v242(
    overlap,
    microstates,
    *,
    minimum_diagonal_overlap=0.5,
):
    """Fix independent real-root signs against the center geometry.

    Polar transport is still used for operator transport.  This narrower operation
    supplies the differentiable adiabatic-root gauge required by the overlap central
    difference used to estimate the derivative connection.  One phase is shared by
    every M_S component belonging to the same spin-free root.
    """

    matrix = np.asarray(overlap, dtype=complex)
    microstates = tuple(microstates)
    if matrix.shape != (len(microstates), len(microstates)):
        raise ValueError("phase alignment requires a square complete-multiplet overlap.")
    root_representatives = {}
    for index, state in enumerate(microstates):
        root_representatives.setdefault(int(state.root_index), index)
    root_phases = {}
    for root_index, representative in root_representatives.items():
        diagonal = matrix[representative, representative]
        magnitude = float(abs(diagonal))
        if magnitude < float(minimum_diagonal_overlap):
            raise ValueError(
                "center-referenced root phase is ambiguous; diagonal overlap is too small."
            )
        root_phases[root_index] = diagonal.conjugate() / magnitude
    phases = np.asarray(
        [root_phases[int(state.root_index)] for state in microstates], dtype=complex
    )
    aligned = matrix @ np.diag(phases)
    return aligned, phases


@dataclass(frozen=True)
class TransportedSOCDerivativeV242:
    coordinate_label: str
    displacement_bohr: float
    center_fingerprint: str
    minus_fingerprint: str
    plus_fingerprint: str
    overlap_center_minus: np.ndarray
    overlap_center_plus: np.ndarray
    transport_minus_to_center: np.ndarray
    transport_plus_to_center: np.ndarray
    H_spin_free_minus_to_center: np.ndarray
    H_spin_free_plus_to_center: np.ndarray
    H_soc_minus_to_center: np.ndarray
    H_soc_plus_to_center: np.ndarray
    K_spin_free: np.ndarray
    K_soc: np.ndarray
    K_total: np.ndarray
    derivative_connection: np.ndarray
    overlap_metrics: dict
    residuals: dict

    def __post_init__(self):
        for name in (
            "overlap_center_minus",
            "overlap_center_plus",
            "transport_minus_to_center",
            "transport_plus_to_center",
            "H_spin_free_minus_to_center",
            "H_spin_free_plus_to_center",
            "H_soc_minus_to_center",
            "H_soc_plus_to_center",
            "K_spin_free",
            "K_soc",
            "K_total",
            "derivative_connection",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=complex).copy())
        object.__setattr__(self, "overlap_metrics", dict(self.overlap_metrics))
        object.__setattr__(self, "residuals", dict(self.residuals))

    def validate(self, tolerance=1.0e-8):
        h = float(self.displacement_bohr)
        if not np.isfinite(h) or h <= 0.0:
            raise ValueError("SOC differential displacement must be finite and positive.")
        if not str(self.coordinate_label).strip():
            raise ValueError("SOC differential coordinate label cannot be empty.")
        fingerprints = (
            self.center_fingerprint,
            self.minus_fingerprint,
            self.plus_fingerprint,
        )
        if any(
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
            for value in fingerprints
        ):
            raise ValueError("SOC differential snapshot fingerprints must be SHA-256 digests.")
        if len(set(fingerprints)) != 3:
            raise ValueError("center, minus, and plus snapshot fingerprints must differ.")
        matrices = (
            self.overlap_center_minus,
            self.overlap_center_plus,
            self.transport_minus_to_center,
            self.transport_plus_to_center,
            self.H_spin_free_minus_to_center,
            self.H_spin_free_plus_to_center,
            self.H_soc_minus_to_center,
            self.H_soc_plus_to_center,
            self.K_spin_free,
            self.K_soc,
            self.K_total,
            self.derivative_connection,
        )
        shape = np.asarray(matrices[0]).shape
        if len(shape) != 2 or shape[0] < 1 or shape[0] != shape[1]:
            raise ValueError("SOC differential matrices must be nonempty and square.")
        if any(np.asarray(matrix).shape != shape for matrix in matrices):
            raise ValueError("SOC differential matrices have incompatible dimensions.")
        if any(not np.all(np.isfinite(matrix)) for matrix in matrices):
            raise ValueError("SOC differential contains non-finite matrix data.")
        for name, matrix in (
            ("K_spin_free", self.K_spin_free),
            ("K_soc", self.K_soc),
            ("K_total", self.K_total),
        ):
            if _scaled_frobenius_v242(matrix, matrix.conj().T) > tolerance:
                raise ValueError(f"{name} is not Hermitian.")
        if _scaled_frobenius_v242(self.K_total, self.K_spin_free + self.K_soc) > tolerance:
            raise ValueError("transported SOC derivative decomposition is inconsistent.")
        expected_spin_free = (
            self.H_spin_free_plus_to_center - self.H_spin_free_minus_to_center
        ) / (2.0 * h)
        expected_soc = (
            self.H_soc_plus_to_center - self.H_soc_minus_to_center
        ) / (2.0 * h)
        if _scaled_frobenius_v242(self.K_spin_free, expected_spin_free) > tolerance:
            raise ValueError("stored spin-free derivative disagrees with transported endpoints.")
        if _scaled_frobenius_v242(self.K_soc, expected_soc) > tolerance:
            raise ValueError("stored SOC derivative disagrees with transported endpoints.")
        if _scaled_frobenius_v242(
            self.derivative_connection, -self.derivative_connection.conj().T
        ) > 5.0e-5:
            raise ValueError("finite-difference derivative connection is not anti-Hermitian.")
        _canonical_v242(self.overlap_metrics)
        _canonical_v242(self.residuals)
        return self

    def as_dict(self, *, include_matrices=True):
        self.validate()
        payload = {
            "coordinate_label": self.coordinate_label,
            "displacement_bohr": float(self.displacement_bohr),
            "center_fingerprint": self.center_fingerprint,
            "minus_fingerprint": self.minus_fingerprint,
            "plus_fingerprint": self.plus_fingerprint,
            "overlap_metrics": self.overlap_metrics,
            "residuals": self.residuals,
        }
        if include_matrices:
            for name in (
                "overlap_center_minus",
                "overlap_center_plus",
                "transport_minus_to_center",
                "transport_plus_to_center",
                "H_spin_free_minus_to_center",
                "H_spin_free_plus_to_center",
                "H_soc_minus_to_center",
                "H_soc_plus_to_center",
                "K_spin_free",
                "K_soc",
                "K_total",
                "derivative_connection",
            ):
                payload[name] = _complex_pairs_v242(getattr(self, name))
        return payload


def transported_soc_central_difference_v242(
    center,
    minus,
    plus,
    *,
    displacement_bohr,
    coordinate_label,
    overlap_policy=FiniteManifoldOverlapPolicyV233(
        minimum_retained_singular_value=0.9,
        maximum_condition_number=10.0,
        maximum_principal_angle_radians=math.acos(0.9),
    ),
):
    """Transport neighboring component matrices to the center before differencing."""

    center = center.validate()
    minus = minus.validate()
    plus = plus.validate()
    h = float(displacement_bohr)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("central-difference displacement must be finite and positive.")
    overlap_minus = complete_multiplet_overlap_v242(center, minus)
    overlap_plus = complete_multiplet_overlap_v242(center, plus)
    transport_minus = certified_transport_from_overlap_v233(
        overlap_minus, policy=overlap_policy
    )
    transport_plus = certified_transport_from_overlap_v233(
        overlap_plus, policy=overlap_policy
    )
    W_minus = transport_minus.right_to_left_transport
    W_plus = transport_plus.right_to_left_transport

    def transported_pair(name):
        left = W_minus @ np.asarray(getattr(minus.matrices, name)) @ W_minus.conj().T
        right = W_plus @ np.asarray(getattr(plus.matrices, name)) @ W_plus.conj().T
        return left, right

    H_spin_free_minus, H_spin_free_plus = transported_pair("H_spin_free")
    H_soc_minus, H_soc_plus = transported_pair("H_soc")
    K_spin_free = (H_spin_free_plus - H_spin_free_minus) / (2.0 * h)
    K_soc = (H_soc_plus - H_soc_minus) / (2.0 * h)
    K_total = (
        (H_spin_free_plus + H_soc_plus)
        - (H_spin_free_minus + H_soc_minus)
    ) / (2.0 * h)
    # The first two OH roots form a degenerate subspace and may undergo an arbitrary
    # orthogonal rotation at each independent CASSCF solve.  Align the complete
    # endpoint spaces with the polar gauge rather than assigning root phases.  The
    # anti-Hermitian part of the aligned overlap slope is the connection in this
    # local parallel-transport gauge; its Hermitian part measures finite-manifold
    # contraction/curvature and is retained as a separate diagnostic.
    aligned_minus = overlap_minus @ W_minus.conj().T
    aligned_plus = overlap_plus @ W_plus.conj().T
    aligned_overlap_slope = (aligned_plus - aligned_minus) / (2.0 * h)
    D = 0.5 * (aligned_overlap_slope - aligned_overlap_slope.conj().T)
    aligned_overlap_hermitian_slope = 0.5 * (
        aligned_overlap_slope + aligned_overlap_slope.conj().T
    )
    identity = np.eye(K_total.shape[0], dtype=complex)
    J = np.asarray(center.matrices.time_reversal_matrix, dtype=complex)
    residuals = {
        "K_spin_free_hermiticity": _scaled_frobenius_v242(
            K_spin_free, K_spin_free.conj().T
        ),
        "K_soc_hermiticity": _scaled_frobenius_v242(K_soc, K_soc.conj().T),
        "K_total_hermiticity": _scaled_frobenius_v242(K_total, K_total.conj().T),
        "K_component_decomposition": _scaled_frobenius_v242(
            K_total, K_spin_free + K_soc
        ),
        "D_antihermiticity": _scaled_frobenius_v242(D, -D.conj().T),
        "minus_transport_unitarity": _scaled_frobenius_v242(
            W_minus.conj().T @ W_minus, identity
        ),
        "plus_transport_unitarity": _scaled_frobenius_v242(
            W_plus.conj().T @ W_plus, identity
        ),
        "K_spin_free_time_reversal": _scaled_frobenius_v242(
            K_spin_free, J @ K_spin_free.conj() @ J.conj().T
        ),
        "K_soc_time_reversal": _scaled_frobenius_v242(
            K_soc, J @ K_soc.conj() @ J.conj().T
        ),
        "K_total_time_reversal": _scaled_frobenius_v242(
            K_total, J @ K_total.conj() @ J.conj().T
        ),
        "parallel_transport_hermitian_slope_frobenius": float(
            np.linalg.norm(aligned_overlap_hermitian_slope, ord="fro")
        ),
    }
    overlap_metrics = {
        "minus": transport_minus.as_dict(),
        "plus": transport_plus.as_dict(),
        "connection_gauge": {
            "convention": (
                "degenerate-safe complete-manifold polar parallel transport; "
                "connection is the anti-Hermitian aligned-overlap slope"
            ),
            "minus_aligned_diagonal_minimum_abs": float(
                np.min(np.abs(np.diag(aligned_minus)))
            ),
            "plus_aligned_diagonal_minimum_abs": float(
                np.min(np.abs(np.diag(aligned_plus)))
            ),
        },
    }
    return TransportedSOCDerivativeV242(
        coordinate_label=str(coordinate_label),
        displacement_bohr=h,
        center_fingerprint=center.fingerprint(),
        minus_fingerprint=minus.fingerprint(),
        plus_fingerprint=plus.fingerprint(),
        overlap_center_minus=overlap_minus,
        overlap_center_plus=overlap_plus,
        transport_minus_to_center=W_minus,
        transport_plus_to_center=W_plus,
        H_spin_free_minus_to_center=H_spin_free_minus,
        H_spin_free_plus_to_center=H_spin_free_plus,
        H_soc_minus_to_center=H_soc_minus,
        H_soc_plus_to_center=H_soc_plus,
        K_spin_free=K_spin_free,
        K_soc=K_soc,
        K_total=K_total,
        derivative_connection=D,
        overlap_metrics=overlap_metrics,
        residuals=residuals,
    ).validate()


@dataclass(frozen=True)
class PySCFSOCDifferentialScanV242:
    runtime_fingerprint: object
    center: PySCFSOCGeometrySnapshotV242
    endpoint_snapshots: tuple
    derivative_records: tuple
    direct_jk_explicit_max_abs_error: float
    convergence_metrics: dict
    claims: dict

    def validate(self):
        runtime = self.runtime_fingerprint.validate()
        center = self.center.validate()
        if runtime.environment_sha256 != center.environment_sha256:
            raise ValueError("runtime and center-snapshot environment identities differ.")
        records = tuple(record.validate() for record in self.derivative_records)
        if len(records) < 3:
            raise ValueError("SOC differential scan requires at least three step sizes.")
        endpoints = tuple(tuple(pair) for pair in self.endpoint_snapshots)
        if len(endpoints) != len(records) or any(len(pair) != 2 for pair in endpoints):
            raise ValueError("one minus/plus snapshot pair is required per derivative record.")
        steps = np.asarray([record.displacement_bohr for record in records], dtype=float)
        if not np.all(steps[:-1] > steps[1:]):
            raise ValueError("SOC differential steps must be strictly decreasing.")
        if any(record.center_fingerprint != center.fingerprint() for record in records):
            raise ValueError("SOC differential records do not share one center snapshot.")
        center_geometry = np.asarray(center.geometry_bohr, dtype=float)
        for record, (minus, plus) in zip(records, endpoints):
            minus = minus.validate()
            plus = plus.validate()
            if (
                minus.environment_sha256 != runtime.environment_sha256
                or plus.environment_sha256 != runtime.environment_sha256
            ):
                raise ValueError("endpoint snapshot runtime identity is inconsistent.")
            if record.minus_fingerprint != minus.fingerprint():
                raise ValueError("minus endpoint fingerprint is not bound to its record.")
            if record.plus_fingerprint != plus.fingerprint():
                raise ValueError("plus endpoint fingerprint is not bound to its record.")
            if minus.state_order != center.state_order or plus.state_order != center.state_order:
                raise ValueError("endpoint complete-multiplet orders differ from the center.")
            if not np.allclose(
                minus.state_average_weights, center.state_average_weights, atol=1.0e-14
            ) or not np.allclose(
                plus.state_average_weights, center.state_average_weights, atol=1.0e-14
            ):
                raise ValueError("endpoint state-average weights differ from the center.")
            expected_minus = center_geometry.copy()
            expected_plus = center_geometry.copy()
            expected_minus[1, 2] -= record.displacement_bohr
            expected_plus[1, 2] += record.displacement_bohr
            if not np.allclose(minus.geometry_bohr, expected_minus, atol=1.0e-12):
                raise ValueError("minus endpoint geometry disagrees with its displacement.")
            if not np.allclose(plus.geometry_bohr, expected_plus, atol=1.0e-12):
                raise ValueError("plus endpoint geometry disagrees with its displacement.")
        error = float(self.direct_jk_explicit_max_abs_error)
        if not np.isfinite(error) or error < 0.0:
            raise ValueError("direct-JK cross-check error must be finite and nonnegative.")
        _canonical_v242(self.convergence_metrics)
        if any(type(value) is not bool for value in self.claims.values()):
            raise TypeError("every v0.24.2 scan claim must be a native Boolean.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242,
            "capability": PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242,
            "runtime": self.runtime_fingerprint.as_dict(),
            "center": self.center.compact_dict(),
            "endpoint_snapshots": [
                {
                    "minus": minus.compact_dict(),
                    "plus": plus.compact_dict(),
                }
                for minus, plus in self.endpoint_snapshots
            ],
            "derivative_records": [record.as_dict() for record in self.derivative_records],
            "direct_jk_explicit_max_abs_error": float(
                self.direct_jk_explicit_max_abs_error
            ),
            "convergence_metrics": self.convergence_metrics,
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v242(self.as_dict())


def _derivative_convergence_metrics_v242(records):
    records = tuple(records)
    metrics = {}
    for name in ("K_spin_free", "K_soc", "K_total", "derivative_connection"):
        changes = [
            float(np.linalg.norm(getattr(fine, name) - getattr(coarse, name), ord="fro"))
            for coarse, fine in zip(records[:-1], records[1:])
        ]
        ratio = float(changes[-1] / changes[-2]) if changes[-2] > 0.0 else 0.0
        richardson = (
            4.0 * np.asarray(getattr(records[-1], name))
            - np.asarray(getattr(records[-2], name))
        ) / 3.0
        metrics[name] = {
            "successive_change_frobenius": changes,
            "fine_to_coarse_change_ratio": ratio,
            "richardson_error_estimate_frobenius": float(changes[-1] / 3.0),
            "finest_norm_frobenius": float(np.linalg.norm(getattr(records[-1], name))),
            "richardson_norm_frobenius": float(np.linalg.norm(richardson)),
        }
    return metrics


def _run_oh_casscf_snapshot_v242(z_hydrogen, *, environment_sha256):
    from pyscf import gto, mcscf, scf

    mol = gto.M(
        atom=(
            ("O", (0.0, 0.0, 0.0)),
            ("H", (0.0, 0.0, float(z_hydrogen))),
        ),
        unit="Bohr",
        basis="sto-3g",
        charge=0,
        spin=1,
        symmetry=False,
        verbose=0,
    )
    mean_field = scf.ROHF(mol)
    mean_field.conv_tol = 1.0e-12
    mean_field.max_cycle = 100
    mean_field.kernel()
    if not mean_field.converged:
        raise RuntimeError("OH ROHF did not converge in the v0.24.2 geometry scan.")
    casscf = mcscf.CASSCF(mean_field, 4, 5).state_average_((1.0 / 3.0,) * 3)
    casscf.conv_tol = 1.0e-9
    casscf.max_cycle_macro = 100
    casscf.kernel()
    if not casscf.converged:
        raise RuntimeError("OH SA-CASSCF did not converge in the v0.24.2 geometry scan.")
    return build_pyscf_soc_geometry_snapshot_v242(
        casscf,
        environment_sha256=environment_sha256,
        root_labels=("D1", "D2", "D3"),
        root_spin_twice=(1, 1, 1),
        weights=(1.0 / 3.0,) * 3,
        molecule_name="OH radical",
        basis_label="STO-3G",
    )


def run_pyscf_oh_bond_differential_soc_v242(
    *,
    steps_bohr=OH_BOND_STEPS_BOHR_V242,
    memory_probe_policy="proc_self",
):
    """Run the connected OH bond-coordinate direct-JK SOC differential scan."""

    steps = tuple(float(step) for step in steps_bohr)
    if len(steps) < 3 or any(not np.isfinite(step) or step <= 0.0 for step in steps):
        raise ValueError("OH SOC scan requires at least three positive finite steps.")
    if any(left <= right for left, right in zip(steps[:-1], steps[1:])):
        raise ValueError("OH SOC scan steps must be strictly decreasing.")
    with guarded_pyscf_runtime_v232(
        memory_probe_policy=memory_probe_policy
    ) as runtime_context:
        environment_sha256 = runtime_context.fingerprint.environment_sha256
        center = _run_oh_casscf_snapshot_v242(
            OH_BOND_LENGTH_BOHR_V242,
            environment_sha256=environment_sha256,
        )
        explicit = build_pyscf_bp_somf_integrals_v241(
            center.wavefunction_snapshot.mol,
            center.wavefunction_snapshot.mo_coeff,
            center.integrals.state_average_density_mo,
        )
        direct_jk_error = float(
            max(
                np.max(
                    np.abs(
                        center.integrals.two_electron_somf_ao_cartesian
                        - explicit.two_electron_somf_ao_cartesian
                    )
                ),
                np.max(
                    np.abs(
                        center.integrals.effective_mo_cartesian
                        - explicit.effective_mo_cartesian
                    )
                ),
            )
        )
        endpoint_snapshots = []
        records = []
        for step in steps:
            minus = _run_oh_casscf_snapshot_v242(
                OH_BOND_LENGTH_BOHR_V242 - step,
                environment_sha256=environment_sha256,
            )
            plus = _run_oh_casscf_snapshot_v242(
                OH_BOND_LENGTH_BOHR_V242 + step,
                environment_sha256=environment_sha256,
            )
            endpoint_snapshots.append((minus, plus))
            records.append(
                transported_soc_central_difference_v242(
                    center,
                    minus,
                    plus,
                    displacement_bohr=step,
                    coordinate_label="H_z",
                )
            )
        metrics = _derivative_convergence_metrics_v242(records)
        claims = {
            "direct_jk_somf_execution_validated": True,
            "rank_five_tensor_avoided_in_production_path": True,
            "connected_geometry_soc_snapshots_validated": True,
            "complete_doublet_overlap_transport_validated": True,
            "degenerate_subspace_polar_gauge_validated": True,
            "transported_spin_free_derivative_preview_validated": True,
            "transported_soc_derivative_preview_validated": True,
            "continuous_physical_derivative_connection_validated": False,
            "full_cartesian_derivative_tensor_validated": False,
            "analytic_soc_derivatives_validated": False,
            "real_mixed_multiplicity_runtime_validated": False,
            "trajectory_ready_molecular_soc_validated": False,
            "live_molecular_soc_backend_admitted": False,
            "ab_initio_soc_accuracy_validated": False,
        }
        return PySCFSOCDifferentialScanV242(
            runtime_fingerprint=runtime_context.fingerprint,
            center=center,
            endpoint_snapshots=tuple(endpoint_snapshots),
            derivative_records=tuple(records),
            direct_jk_explicit_max_abs_error=direct_jk_error,
            convergence_metrics=metrics,
            claims=claims,
        ).validate()


@dataclass(frozen=True)
class PySCFSOCDifferentialAuditV242:
    checks: dict
    metrics: dict
    thresholds: dict
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or not self.checks:
            raise ValueError("v0.24.2 differential audit requires checks.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.24.2 differential audit gate must be Boolean.")
        _canonical_v242(self.metrics)
        _canonical_v242(self.thresholds)
        if type(self.passed) is not bool:
            raise TypeError("v0.24.2 differential audit result must be Boolean.")
        if self.passed != bool(all(self.checks.values())):
            raise ValueError("v0.24.2 differential audit result disagrees with its gates.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "checks": dict(self.checks),
            "metrics": self.metrics,
            "thresholds": self.thresholds,
            "passed": bool(self.passed),
        }


def audit_pyscf_oh_bond_differential_soc_v242(scan):
    """Audit the real direct-JK connected-geometry OH evidence."""

    scan = scan.validate()
    center = scan.center
    records = scan.derivative_records
    endpoints = scan.endpoint_snapshots
    runtime = scan.runtime_fingerprint
    thresholds = {
        "direct_jk_explicit_max_abs_error": 1.0e-12,
        "minimum_overlap_singular_value": 0.99,
        "matrix_residual": 1.0e-9,
        "derivative_connection_residual": 1.0e-12,
        "minimum_second_order_ratio": 0.15,
        "maximum_second_order_ratio": 0.40,
        "maximum_relative_richardson_error": 5.0e-3,
        "minimum_nonzero_soc_derivative_norm": 1.0e-8,
        "maximum_kramers_splitting_hartree": 1.0e-10,
    }
    all_overlap_blocks = [
        record.overlap_metrics[orientation]
        for record in records
        for orientation in ("minus", "plus")
    ]
    residual_names = (
        "K_spin_free_hermiticity",
        "K_soc_hermiticity",
        "K_total_hermiticity",
        "K_component_decomposition",
        "minus_transport_unitarity",
        "plus_transport_unitarity",
        "K_spin_free_time_reversal",
        "K_soc_time_reversal",
        "K_total_time_reversal",
    )
    maximum_matrix_residual = float(
        max(record.residuals[name] for record in records for name in residual_names)
    )
    maximum_D_residual = float(
        max(record.residuals["D_antihermiticity"] for record in records)
    )
    minimum_overlap_singular = float(
        min(block["minimum_singular_value"] for block in all_overlap_blocks)
    )
    maximum_overlap_singular = float(
        max(block["maximum_singular_value"] for block in all_overlap_blocks)
    )
    maximum_transport_residual = float(
        max(block["transport_unitarity_residual"] for block in all_overlap_blocks)
    )
    hermitian_slopes = [
        float(record.residuals["parallel_transport_hermitian_slope_frobenius"])
        for record in records
    ]
    relative_richardson = {}
    for name in ("K_spin_free", "K_soc", "K_total"):
        item = scan.convergence_metrics[name]
        relative_richardson[name] = float(
            item["richardson_error_estimate_frobenius"]
            / max(item["richardson_norm_frobenius"], 1.0e-30)
        )
    ratios = {
        name: float(scan.convergence_metrics[name]["fine_to_coarse_change_ratio"])
        for name in ("K_spin_free", "K_soc", "K_total")
    }
    expected_order = (
        "D1(M=+1/2)",
        "D1(M=-1/2)",
        "D2(M=+1/2)",
        "D2(M=-1/2)",
        "D3(M=+1/2)",
        "D3(M=-1/2)",
    )
    thread_values = tuple(runtime.thread_environment.values())
    checks = {
        "runtime_exact_pyscf_distribution": bool(
            runtime.pyscf_distribution_version == "2.13.1"
        ),
        "runtime_exact_pyscf_module": bool(runtime.pyscf_module_version == "2.13.1"),
        "runtime_numpy_locked": bool(runtime.numpy_version == "2.5.2"),
        "runtime_scipy_locked": bool(runtime.scipy_version == "1.18.0"),
        "runtime_h5py_locked": bool(runtime.h5py_version == "3.16.0"),
        "runtime_threads_fixed_to_one": bool(
            thread_values and all(str(value) == "1" for value in thread_values)
        ),
        "runtime_identity_matches_snapshots": bool(
            runtime.environment_sha256 == center.environment_sha256
        ),
        "all_six_endpoint_receipts_are_retained": bool(
            len(endpoints) == 3 and all(len(pair) == 2 for pair in endpoints)
        ),
        "endpoint_runtime_identities_match": bool(
            all(
                snapshot.environment_sha256 == runtime.environment_sha256
                for pair in endpoints
                for snapshot in pair
            )
        ),
        "endpoint_fingerprints_are_bound_to_records": bool(
            all(
                record.minus_fingerprint == pair[0].fingerprint()
                and record.plus_fingerprint == pair[1].fingerprint()
                for record, pair in zip(records, endpoints)
            )
        ),
        "endpoint_geometries_match_signed_displacements": bool(
            all(
                abs(pair[0].geometry_bohr[1, 2] - (OH_BOND_LENGTH_BOHR_V242 - record.displacement_bohr)) <= 1.0e-12
                and abs(pair[1].geometry_bohr[1, 2] - (OH_BOND_LENGTH_BOHR_V242 + record.displacement_bohr)) <= 1.0e-12
                for record, pair in zip(records, endpoints)
            )
        ),
        "direct_jk_strategy_declared": bool(
            center.calculation_input["somf_contraction"]
            == PYSCF_DIRECT_JK_SOMF_STRATEGY_V242
        ),
        "direct_jk_matches_explicit_tensor_oracle": bool(
            scan.direct_jk_explicit_max_abs_error
            <= thresholds["direct_jk_explicit_max_abs_error"]
        ),
        "one_electron_integral_identity_frozen": bool(
            center.integrals.one_electron_integral
            == BP_SOMF_ONE_ELECTRON_INTEGRAL_V241
        ),
        "two_electron_integral_identity_frozen": bool(
            center.integrals.two_electron_integral
            == BP_SOMF_TWO_ELECTRON_INTEGRAL_V241
        ),
        "single_half_over_c_squared_prefactor": bool(
            abs(
                center.integrals.prefactor
                - 0.5 / center.integrals.light_speed_au**2
            )
            <= 1.0e-16
        ),
        "state_average_density_has_nine_electrons": bool(
            abs(float(np.trace(center.integrals.state_average_density_mo).real) - 9.0)
            <= 1.0e-10
        ),
        "central_molecular_soc_signal_nonzero": bool(
            np.linalg.norm(center.matrices.H_soc, ord="fro") > 1.0e-6
        ),
        "three_spin_free_doublet_roots": bool(
            len(center.roots) == 3
            and all(root.spin_twice == 1 for root in center.roots)
        ),
        "six_complete_doublet_microstates": bool(len(center.state_order) == 6),
        "exact_complete_microstate_order": bool(center.state_order == expected_order),
        "all_center_roots_spin_pure": bool(
            all(
                abs(float(root.spin_square) - 0.75) <= 1.0e-6
                for root in center.roots
            )
        ),
        "center_soc_matrix_hermitian": bool(
            _scaled_frobenius_v242(
                center.matrices.H_soc, center.matrices.H_soc.conj().T
            )
            <= thresholds["matrix_residual"]
        ),
        "center_total_time_reversal_invariant": bool(
            center.matrices.time_reversal_residual <= thresholds["matrix_residual"]
        ),
        "center_time_reversal_square_is_fermionic": bool(
            center.matrices.time_reversal_square_residual
            <= thresholds["matrix_residual"]
        ),
        "center_kramers_pairs_resolved": bool(
            center.matrices.maximum_kramers_pair_splitting_hartree
            <= thresholds["maximum_kramers_splitting_hartree"]
        ),
        "all_raw_overlaps_are_physical_contractions": bool(
            all(block["physically_consistent"] for block in all_overlap_blocks)
        ),
        "all_polar_transports_are_trajectory_quality": bool(
            all(block["trajectory_ready"] for block in all_overlap_blocks)
        ),
        "all_overlap_singular_values_are_retained": bool(
            minimum_overlap_singular >= thresholds["minimum_overlap_singular_value"]
        ),
        "no_overlap_spectral_expansion": bool(maximum_overlap_singular <= 1.0 + 1.0e-10),
        "all_polar_transports_are_unitary": bool(
            maximum_transport_residual <= thresholds["matrix_residual"]
        ),
        "all_polar_factors_are_positive": bool(
            all(block["polar_minimum_eigenvalue"] >= -1.0e-10 for block in all_overlap_blocks)
        ),
        "parallel_transport_contraction_slope_refines": bool(
            hermitian_slopes[-1] < hermitian_slopes[0]
        ),
        "three_strictly_decreasing_displacements": bool(len(records) == 3),
        "displacements_halve_exactly": bool(
            all(
                abs(records[index + 1].displacement_bohr / records[index].displacement_bohr - 0.5)
                <= 1.0e-14
                for index in range(len(records) - 1)
            )
        ),
        "canonical_displacement_ladder_used": bool(
            tuple(record.displacement_bohr for record in records)
            == OH_BOND_STEPS_BOHR_V242
        ),
        "all_spin_free_derivatives_hermitian": bool(
            all(record.residuals["K_spin_free_hermiticity"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_soc_derivatives_hermitian": bool(
            all(record.residuals["K_soc_hermiticity"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_total_derivatives_hermitian": bool(
            all(record.residuals["K_total_hermiticity"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_derivative_component_sums_exact": bool(
            all(record.residuals["K_component_decomposition"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_spin_free_derivatives_time_reversal_invariant": bool(
            all(record.residuals["K_spin_free_time_reversal"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_soc_derivatives_time_reversal_invariant": bool(
            all(record.residuals["K_soc_time_reversal"] <= thresholds["matrix_residual"] for record in records)
        ),
        "all_total_derivatives_time_reversal_invariant": bool(
            all(record.residuals["K_total_time_reversal"] <= thresholds["matrix_residual"] for record in records)
        ),
        "polar_gauge_connections_antihermitian": bool(
            maximum_D_residual <= thresholds["derivative_connection_residual"]
        ),
        "spin_free_difference_on_second_order_plateau": bool(
            thresholds["minimum_second_order_ratio"]
            <= ratios["K_spin_free"]
            <= thresholds["maximum_second_order_ratio"]
        ),
        "soc_difference_on_second_order_plateau": bool(
            thresholds["minimum_second_order_ratio"]
            <= ratios["K_soc"]
            <= thresholds["maximum_second_order_ratio"]
        ),
        "total_difference_on_second_order_plateau": bool(
            thresholds["minimum_second_order_ratio"]
            <= ratios["K_total"]
            <= thresholds["maximum_second_order_ratio"]
        ),
        "spin_free_richardson_error_bounded": bool(
            relative_richardson["K_spin_free"]
            <= thresholds["maximum_relative_richardson_error"]
        ),
        "soc_richardson_error_bounded": bool(
            relative_richardson["K_soc"]
            <= thresholds["maximum_relative_richardson_error"]
        ),
        "total_richardson_error_bounded": bool(
            relative_richardson["K_total"]
            <= thresholds["maximum_relative_richardson_error"]
        ),
        "finest_soc_derivative_is_nonzero": bool(
            scan.convergence_metrics["K_soc"]["finest_norm_frobenius"]
            >= thresholds["minimum_nonzero_soc_derivative_norm"]
        ),
        "maximum_matrix_residual_bounded": bool(
            maximum_matrix_residual <= thresholds["matrix_residual"]
        ),
        "capability_is_differential_preview": bool(
            PYSCF_DIFFERENTIAL_SOC_CAPABILITY_V242
            == "connected_geometry_differential_preview"
        ),
        "full_cartesian_derivative_claim_remains_false": bool(
            scan.claims["full_cartesian_derivative_tensor_validated"] is False
        ),
        "analytic_soc_derivative_claim_remains_false": bool(
            scan.claims["analytic_soc_derivatives_validated"] is False
        ),
        "continuous_connection_claim_remains_false": bool(
            scan.claims["continuous_physical_derivative_connection_validated"] is False
        ),
        "real_mixed_multiplicity_claim_remains_false": bool(
            scan.claims["real_mixed_multiplicity_runtime_validated"] is False
        ),
        "trajectory_ready_claim_remains_false": bool(
            scan.claims["trajectory_ready_molecular_soc_validated"] is False
        ),
        "live_backend_admission_remains_false": bool(
            scan.claims["live_molecular_soc_backend_admitted"] is False
        ),
        "accuracy_claim_remains_false": bool(
            scan.claims["ab_initio_soc_accuracy_validated"] is False
        ),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    metrics = {
        "direct_jk_explicit_max_abs_error": float(
            scan.direct_jk_explicit_max_abs_error
        ),
        "minimum_overlap_singular_value": minimum_overlap_singular,
        "maximum_overlap_singular_value": maximum_overlap_singular,
        "maximum_transport_unitarity_residual": maximum_transport_residual,
        "maximum_matrix_residual": maximum_matrix_residual,
        "maximum_D_antihermiticity_residual": maximum_D_residual,
        "parallel_transport_hermitian_slopes": hermitian_slopes,
        "second_order_ratios": ratios,
        "relative_richardson_errors": relative_richardson,
        "finest_K_spin_free_norm": float(
            scan.convergence_metrics["K_spin_free"]["finest_norm_frobenius"]
        ),
        "finest_K_soc_norm": float(
            scan.convergence_metrics["K_soc"]["finest_norm_frobenius"]
        ),
        "finest_K_total_norm": float(
            scan.convergence_metrics["K_total"]["finest_norm_frobenius"]
        ),
        "gate_count": len(checks),
        "retained_endpoint_snapshot_count": int(
            sum(len(pair) for pair in endpoints)
        ),
    }
    return PySCFSOCDifferentialAuditV242(
        checks=checks,
        metrics=metrics,
        thresholds=thresholds,
        passed=bool(all(checks.values())),
    ).validate()


@dataclass(frozen=True)
class PySCFSOCDifferentialEvidenceV242:
    scan: PySCFSOCDifferentialScanV242
    audit: PySCFSOCDifferentialAuditV242

    def validate(self):
        self.scan.validate()
        self.audit.validate()
        if not self.audit.passed:
            failed = ", ".join(
                name for name, passed in self.audit.checks.items() if not passed
            )
            raise ValueError("v0.24.2 differential evidence failed: " + failed)
        return self

    @property
    def claims(self):
        return dict(self.scan.claims)

    def as_dict(self):
        self.validate()
        return {
            "schema": PYSCF_DIFFERENTIAL_SOC_SCHEMA_V242,
            "scan": self.scan.as_dict(),
            "audit": self.audit.as_dict(),
            "claims": self.claims,
        }

    def fingerprint(self):
        return _sha256_v242(self.as_dict())


def run_pyscf_oh_bond_differential_evidence_v242(**kwargs):
    scan = run_pyscf_oh_bond_differential_soc_v242(**kwargs)
    return PySCFSOCDifferentialEvidenceV242(
        scan=scan,
        audit=audit_pyscf_oh_bond_differential_soc_v242(scan),
    ).validate()


def save_pyscf_oh_bond_differential_evidence_v242(path, evidence=None):
    from .campaign_io import save_campaign_json

    evidence = (
        run_pyscf_oh_bond_differential_evidence_v242()
        if evidence is None
        else evidence.validate()
    )
    return save_campaign_json(path, evidence.as_dict())
