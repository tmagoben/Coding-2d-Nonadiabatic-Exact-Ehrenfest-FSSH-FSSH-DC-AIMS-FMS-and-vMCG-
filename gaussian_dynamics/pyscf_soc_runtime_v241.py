"""Pinned real-runtime evidence for the v0.24.1 PySCF static SOC provider."""

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from .electronic_contract_v213 import hartree_to_wavenumber_v213
from .pyscf_runtime_v232 import guarded_pyscf_runtime_v232
from .pyscf_state_interaction_soc_v241 import (
    BP_SOMF_ONE_ELECTRON_INTEGRAL_V241,
    BP_SOMF_OPERATOR_FAMILY_V241,
    BP_SOMF_STATIC_LIMITATION_V241,
    BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    PYSCF_REQUIRED_VERSION_V241,
    PySCFStateInteractionSOCProviderV241,
    PySCFStateInteractionSOCResultV241,
    complete_spin_microstates_v241,
)


PYSCF_SOC_RUNTIME_SCHEMA_V241 = "gnd-pyscf-static-soc-runtime-evidence-v0.24.1"
OH_BOND_LENGTH_BOHR_V241 = 1.83256418024373
OH_ISOTOPE_MASSES_AMU_V241 = (15.99491461957, 1.00782503223)


def _canonical_runtime_v241(value):
    if isinstance(value, np.generic):
        return _canonical_runtime_v241(value.item())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("runtime evidence cannot contain non-finite data.")
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.ndarray):
        return _canonical_runtime_v241(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("runtime-evidence dictionary keys must be strings.")
        return {
            key: _canonical_runtime_v241(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_runtime_v241(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("runtime evidence cannot contain non-finite data.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported runtime-evidence value {type(value).__name__}.")


def _canonical_runtime_bytes_v241(value):
    return json.dumps(
        _canonical_runtime_v241(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _scaled_error_v241(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def crosscheck_pyscf_somf_jk_v241(mol, density_ao, expected_two_electron_somf):
    """Cross-check the explicit rank-five contraction with PySCF's JK driver."""

    from pyscf.scf import jk

    density = np.asarray(density_ao, dtype=float)
    expected = np.asarray(expected_two_electron_somf, dtype=float)
    direct, exchange_left, exchange_right = jk.get_jk(
        mol,
        [density, density, density],
        scripts=(
            "ijkl,kl->ij",
            "ijkl,jk->il",
            "ijkl,li->kj",
        ),
        intor=BP_SOMF_TWO_ELECTRON_INTEGRAL_V241,
    )
    observed = (
        np.asarray(direct)
        - 1.5 * np.asarray(exchange_left)
        - 1.5 * np.asarray(exchange_right)
    )
    if observed.shape != expected.shape or not np.all(np.isfinite(observed)):
        raise ValueError("PySCF JK SOC cross-check returned invalid data.")
    return float(np.max(np.abs(observed - expected)))


@dataclass(frozen=True)
class PySCFStaticSOCAuditV241:
    checks: dict
    metrics: dict
    static_soc_validated: bool
    trajectory_ready: bool
    live_backend_admitted: bool
    ab_initio_accuracy_validated: bool
    passed: bool

    def validate(self):
        if not isinstance(self.checks, dict) or not self.checks:
            raise ValueError("static SOC audit requires checks.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every static SOC audit gate must be a native Boolean.")
        if not isinstance(self.metrics, dict):
            raise TypeError("static SOC audit metrics must be a dictionary.")
        _canonical_runtime_v241(self.metrics)
        for name, value in (
            ("static_soc_validated", self.static_soc_validated),
            ("trajectory_ready", self.trajectory_ready),
            ("live_backend_admitted", self.live_backend_admitted),
            ("ab_initio_accuracy_validated", self.ab_initio_accuracy_validated),
            ("passed", self.passed),
        ):
            if type(value) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if self.passed != all(self.checks.values()):
            raise ValueError("static SOC audit pass flag disagrees with its gates.")
        if not self.static_soc_validated or not self.passed:
            raise ValueError("recorded v0.24.1 static SOC audit did not pass.")
        if self.trajectory_ready or self.live_backend_admitted:
            raise ValueError("static-only PySCF SOC cannot be trajectory admitted.")
        if self.ab_initio_accuracy_validated:
            raise ValueError(
                "one STO-3G smoke calculation cannot establish ab-initio accuracy."
            )
        return self

    def as_dict(self):
        self.validate()
        return {
            "checks": dict(self.checks),
            "metrics": _canonical_runtime_v241(self.metrics),
            "static_soc_validated": self.static_soc_validated,
            "trajectory_ready": self.trajectory_ready,
            "live_backend_admitted": self.live_backend_admitted,
            "ab_initio_accuracy_validated": self.ab_initio_accuracy_validated,
            "passed": self.passed,
        }


def audit_pyscf_static_soc_v241(result, *, somf_jk_crosscheck_error):
    """Audit a real fixed-geometry result without promoting dynamics claims."""

    if not isinstance(result, PySCFStateInteractionSOCResultV241):
        raise TypeError("static SOC audit requires PySCFStateInteractionSOCResultV241.")
    result = result.validate()
    matrices = result.matrices
    integrals = result.integrals
    roots = matrices.roots
    microstates = matrices.microstates
    H0 = matrices.H_spin_free
    Hso = matrices.H_soc
    H = matrices.H_total
    eigenvectors = matrices.soc_eigenvectors
    eigenvalues = matrices.soc_eigenvalues_hartree
    nstate = len(microstates)
    identity = np.eye(nstate, dtype=complex)
    electron_count = int(result.identity.electron_count)
    density_trace_error = abs(
        float(np.trace(integrals.state_average_density_mo).real) - electron_count
    )
    eigensystem_error = _scaled_error_v241(
        H,
        eigenvectors @ np.diag(eigenvalues) @ eigenvectors.conj().T,
    )
    component_error = _scaled_error_v241(H, H0 + Hso)
    H0_hermiticity = _scaled_error_v241(H0, H0.conj().T)
    Hso_hermiticity = _scaled_error_v241(Hso, Hso.conj().T)
    H_hermiticity = _scaled_error_v241(H, H.conj().T)
    eigenvector_unitarity = _scaled_error_v241(
        eigenvectors.conj().T @ eigenvectors, identity
    )
    complete_count = sum(root.multiplicity for root in roots)
    root_spin_errors = [
        abs(
            float(root.spin_square)
            - 0.25 * root.spin_twice * (root.spin_twice + 2)
        )
        for root in roots
    ]
    dynamic_flags_false = bool(
        not result.capabilities.spin_free_derivatives
        and not result.capabilities.soc_derivatives
        and not result.capabilities.derivative_connections
        and not result.capabilities.cross_geometry_overlaps
        and not result.capabilities.analytic_soc_derivatives
    )
    odd_electron = result.identity.electron_parity == "odd"
    if odd_electron:
        kramers_gate = bool(
            matrices.maximum_kramers_pair_splitting_hartree is not None
            and matrices.maximum_kramers_pair_splitting_hartree <= 1.0e-10
        )
    else:
        kramers_gate = matrices.maximum_kramers_pair_splitting_hartree is None

    checks = {
        "exact_pyscf_distribution_and_module_version": bool(
            result.runtime_probe.exact_version
            and result.runtime_probe.distribution_version
            == PYSCF_REQUIRED_VERSION_V241
            and result.runtime_probe.module_version == PYSCF_REQUIRED_VERSION_V241
        ),
        "required_bp_integral_apis_present": bool(
            result.runtime_probe.integral_apis_available
        ),
        "spin_separated_transition_rdm_api_present": bool(
            result.runtime_probe.transition_rdm_api_available
        ),
        "determinant_spin_ladder_apis_present": bool(
            result.runtime_probe.spin_ladder_apis_available
        ),
        "scf_converged": result.scf_converged is True,
        "casscf_converged": result.casscf_converged is True,
        "direct_soc_assembly_completed": result.soc_assembled is True,
        "live_ab_initio_source_declared": result.identity.source_kind
        == "live_ab_initio",
        "traceable_nuclear_and_runtime_identity": bool(
            result.identity.traceable_nuclear_identity
        ),
        "operator_is_exactly_bp_somf": result.identity.soc_operator
        == BP_SOMF_OPERATOR_FAMILY_V241,
        "one_electron_integral_identity_frozen": (
            integrals.one_electron_integral
            == BP_SOMF_ONE_ELECTRON_INTEGRAL_V241
        ),
        "two_electron_integral_identity_frozen": (
            integrals.two_electron_integral
            == BP_SOMF_TWO_ELECTRON_INTEGRAL_V241
        ),
        "single_prefactor_is_half_over_c_squared": abs(
            integrals.prefactor - 0.5 / integrals.light_speed_au**2
        )
        <= 1.0e-16,
        "state_average_density_has_correct_electron_count": density_trace_error
        <= 1.0e-10,
        "one_electron_ao_soc_is_antisymmetric": (
            integrals.one_electron_antisymmetry_residual <= 1.0e-10
        ),
        "two_electron_somf_ao_is_antisymmetric": (
            integrals.two_electron_antisymmetry_residual <= 1.0e-10
        ),
        "explicit_somf_matches_independent_pyscf_jk_path": float(
            somf_jk_crosscheck_error
        )
        <= 1.0e-12,
        "complete_microstate_count": complete_count == nstate,
        "complete_microstate_order_reconstructs": tuple(microstates)
        == complete_spin_microstates_v241(roots),
        "microstate_labels_are_unique": len(set(matrices.state_order)) == nstate,
        "all_roots_are_spin_pure": max(root_spin_errors, default=0.0) <= 1.0e-6,
        "spin_free_matrix_is_hermitian": H0_hermiticity <= 1.0e-12,
        "direct_soc_matrix_is_hermitian": Hso_hermiticity <= 1.0e-12,
        "total_matrix_is_hermitian": H_hermiticity <= 1.0e-12,
        "total_matrix_decomposition_is_exact": component_error <= 1.0e-12,
        "soc_eigenvectors_are_unitary": eigenvector_unitarity <= 1.0e-12,
        "soc_eigensystem_reconstructs_total_matrix": eigensystem_error <= 1.0e-12,
        "molecular_soc_signal_is_nonzero": float(np.linalg.norm(Hso)) > 1.0e-12,
        "time_reversal_invariance": matrices.time_reversal_residual <= 1.0e-12,
        "time_reversal_square_matches_electron_parity": (
            matrices.time_reversal_square_residual <= 1.0e-12
        ),
        "kramers_or_even_sector_condition": kramers_gate,
        "capability_tier_is_static_soc": result.capabilities.tier == "static_soc",
        "all_unimplemented_dynamic_capabilities_are_false": dynamic_flags_false,
        "trajectory_ready_is_false": result.trajectory_ready is False,
        "v230_real_backend_admission_remains_closed": (
            not result.molecular_soc_contract.real_backend_admission_ready
        ),
        "convention_state_order_matches_direct_matrix": tuple(
            result.convention.state_order
        )
        == matrices.state_order,
        "convention_excludes_external_magnetic_field": (
            result.convention.external_magnetic_field is False
        ),
        "static_limitation_is_explicit_in_identity": (
            result.identity.derivative_method == BP_SOMF_STATIC_LIMITATION_V241
        ),
        "nac_convention_not_reinterpreted_by_static_soc": (
            "not exercised by static SOC"
            in result.provenance.parameters["nac_convention"]
        ),
    }
    checks = {name: bool(value) for name, value in checks.items()}
    metrics = {
        "somf_jk_crosscheck_max_abs_error": float(somf_jk_crosscheck_error),
        "state_average_density_trace_error": float(density_trace_error),
        "H_spin_free_hermiticity_scaled_error": H0_hermiticity,
        "H_soc_hermiticity_scaled_error": Hso_hermiticity,
        "H_total_hermiticity_scaled_error": H_hermiticity,
        "H_total_component_scaled_error": component_error,
        "eigensystem_reconstruction_scaled_error": eigensystem_error,
        "eigenvector_unitarity_scaled_error": eigenvector_unitarity,
        "H_soc_frobenius_norm_hartree": float(np.linalg.norm(Hso)),
        "H_soc_frobenius_norm_cm_inverse": float(
            hartree_to_wavenumber_v213(np.linalg.norm(Hso))
        ),
        "H_soc_max_abs_hartree": float(np.max(np.abs(Hso))),
        "H_soc_max_abs_cm_inverse": float(
            hartree_to_wavenumber_v213(np.max(np.abs(Hso)))
        ),
        "H_soc_trace_abs_hartree": float(abs(np.trace(Hso))),
        "maximum_root_spin_square_error": float(max(root_spin_errors, default=0.0)),
        "time_reversal_residual": float(matrices.time_reversal_residual),
        "time_reversal_square_residual": float(
            matrices.time_reversal_square_residual
        ),
        "maximum_kramers_pair_splitting_hartree": (
            matrices.maximum_kramers_pair_splitting_hartree
        ),
        "maximum_kramers_pair_splitting_cm_inverse": (
            None
            if matrices.maximum_kramers_pair_splitting_hartree is None
            else float(
                hartree_to_wavenumber_v213(
                    matrices.maximum_kramers_pair_splitting_hartree
                )
            )
        ),
        "nroot": len(roots),
        "n_microstate": nstate,
        "nmo": integrals.effective_mo_cartesian.shape[1],
    }
    return PySCFStaticSOCAuditV241(
        checks=checks,
        metrics=metrics,
        static_soc_validated=bool(all(checks.values())),
        trajectory_ready=False,
        live_backend_admitted=False,
        ab_initio_accuracy_validated=False,
        passed=bool(all(checks.values())),
    ).validate()


@dataclass(frozen=True)
class PySCFStaticSOCRuntimeEvidenceV241:
    schema: str
    runtime: object
    calculation: dict
    result: PySCFStateInteractionSOCResultV241
    audit: PySCFStaticSOCAuditV241
    claims: dict

    def validate(self):
        if self.schema != PYSCF_SOC_RUNTIME_SCHEMA_V241:
            raise ValueError("PySCF static SOC runtime schema mismatch.")
        self.runtime.validate()
        self.result.validate()
        self.audit.validate()
        if not isinstance(self.calculation, dict):
            raise TypeError("runtime calculation specification must be a dictionary.")
        _canonical_runtime_v241(self.calculation)
        if not isinstance(self.claims, dict) or any(
            type(value) is not bool for value in self.claims.values()
        ):
            raise TypeError("runtime claims must be a dictionary of native Booleans.")
        required_claims = {
            "real_PySCF_BP_SOMF_execution_validated": True,
            "direct_molecular_SOC_elements_returned": True,
            "doublet_and_Kramers_sector_validated": True,
            "static_molecular_SOC_tier_validated": True,
            "trajectory_ready_molecular_SOC_validated": False,
            "live_molecular_SOC_backend_admitted": False,
            "physical_SOC_derivatives_validated": False,
            "cross_geometry_SOC_tracking_validated": False,
            "ab_initio_SOC_accuracy_validated": False,
            "Prism_runtime_dependency_required": False,
        }
        if self.claims != required_claims:
            raise ValueError("PySCF static SOC runtime claim boundary changed.")
        return self

    @property
    def passed(self):
        return self.audit.passed

    def as_dict(self):
        self.validate()
        return {
            "schema": self.schema,
            "runtime": self.runtime.as_dict(),
            "calculation": _canonical_runtime_v241(self.calculation),
            "result": self.result.as_dict(include_large_arrays=False),
            "audit": self.audit.as_dict(),
            "claims": dict(self.claims),
        }

    def fingerprint(self):
        return hashlib.sha256(
            _canonical_runtime_bytes_v241(self.as_dict())
        ).hexdigest()


def run_pyscf_oh_static_soc_evidence_v241(*, memory_probe_policy="proc_self"):
    """Run the pinned OH SA-CASSCF(5e,4o)/STO-3G doublet evidence case."""

    with guarded_pyscf_runtime_v232(
        memory_probe_policy=memory_probe_policy
    ) as runtime_context:
        from pyscf import gto, mcscf, scf

        calculation = {
            "molecule": "OH radical",
            "atom_symbols": ["O", "H"],
            "geometry_bohr": [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, OH_BOND_LENGTH_BOHR_V241],
            ],
            "isotope_masses_amu": list(OH_ISOTOPE_MASSES_AMU_V241),
            "charge": 0,
            "spin_twice": 1,
            "basis": "STO-3G",
            "scf_method": "ROHF",
            "casscf_method": "equal-weight three-root SA-CASSCF(5e,4o)",
            "root_labels": ["D1", "D2", "D3"],
            "root_spin_twice": [1, 1, 1],
            "state_average_weights": [1.0 / 3.0] * 3,
            "soc_operator": BP_SOMF_OPERATOR_FAMILY_V241,
            "scalar_relativistic_method": "none",
            "scf_convergence_tolerance": 1.0e-12,
            "casscf_convergence_tolerance": 1.0e-9,
            "thread_count": 1,
            "purpose": (
                "fixed-geometry implementation evidence, not a basis/method "
                "convergence or spectroscopic-accuracy benchmark"
            ),
        }
        mol = gto.M(
            atom=(
                ("O", (0.0, 0.0, 0.0)),
                ("H", (0.0, 0.0, OH_BOND_LENGTH_BOHR_V241)),
            ),
            unit="Bohr",
            basis="sto-3g",
            charge=0,
            spin=1,
            symmetry=False,
            verbose=0,
        )
        mean_field = scf.ROHF(mol)
        mean_field.conv_tol = calculation["scf_convergence_tolerance"]
        mean_field.max_cycle = 100
        mean_field.kernel()
        if not mean_field.converged:
            raise RuntimeError("OH ROHF did not converge.")

        casscf = mcscf.CASSCF(mean_field, 4, 5).state_average_(
            calculation["state_average_weights"]
        )
        casscf.conv_tol = calculation["casscf_convergence_tolerance"]
        casscf.max_cycle_macro = 100
        casscf.kernel()
        if not casscf.converged:
            raise RuntimeError("OH three-root SA-CASSCF did not converge.")

        provider = PySCFStateInteractionSOCProviderV241(
            casscf,
            environment_sha256=runtime_context.fingerprint.environment_sha256,
            root_labels=tuple(calculation["root_labels"]),
            root_spin_twice=tuple(calculation["root_spin_twice"]),
            weights=tuple(calculation["state_average_weights"]),
            molecule_name=calculation["molecule"],
            basis_label=calculation["basis"],
            isotope_masses_amu=OH_ISOTOPE_MASSES_AMU_V241,
        )
        result = provider.evaluate_static_soc()
        somf_jk_error = crosscheck_pyscf_somf_jk_v241(
            mol,
            result.integrals.state_average_density_ao,
            result.integrals.two_electron_somf_ao_cartesian,
        )
        audit = audit_pyscf_static_soc_v241(
            result, somf_jk_crosscheck_error=somf_jk_error
        )
        claims = {
            "real_PySCF_BP_SOMF_execution_validated": True,
            "direct_molecular_SOC_elements_returned": True,
            "doublet_and_Kramers_sector_validated": True,
            "static_molecular_SOC_tier_validated": True,
            "trajectory_ready_molecular_SOC_validated": False,
            "live_molecular_SOC_backend_admitted": False,
            "physical_SOC_derivatives_validated": False,
            "cross_geometry_SOC_tracking_validated": False,
            "ab_initio_SOC_accuracy_validated": False,
            "Prism_runtime_dependency_required": False,
        }
        return PySCFStaticSOCRuntimeEvidenceV241(
            schema=PYSCF_SOC_RUNTIME_SCHEMA_V241,
            runtime=runtime_context.fingerprint,
            calculation=calculation,
            result=result,
            audit=audit,
            claims=claims,
        ).validate()


def save_pyscf_oh_static_soc_evidence_v241(path, evidence=None):
    evidence = (
        run_pyscf_oh_static_soc_evidence_v241()
        if evidence is None
        else evidence.validate()
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        evidence.as_dict(),
        sort_keys=True,
        indent=2,
        allow_nan=False,
    )
    path.write_text(payload + "\n", encoding="utf-8")
    return path
