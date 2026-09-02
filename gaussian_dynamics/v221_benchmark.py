"""Corrective hardening and backend-readiness campaign for v0.22.1."""

import copy
from dataclasses import asdict, dataclass, replace
import importlib.util
import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    SOCOperatorComponentsV220,
)
from .block_dynamics_v21 import (
    PrescribedBlockDynamicsSettingsV21,
    run_prescribed_block_dynamics_v21,
)
from .block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,
    BlockSparseSettingsV21,
)
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from .electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
    compose_electronic_operator_v213,
)
from .electronic_operator_v21 import ElectronicOperatorSnapshotV21
from .initial_projection_v213 import project_grid_wavefunction_fixed_frame_v213
from .physical_soc_validation_v220 import audit_physical_soc_provider_v220
from .soc_admission_v221 import (
    SOCSymmetryContractV221,
    audit_soc_symmetry_contract_v221,
)
from .spinor_exact_grid_v220 import (
    SpinorGridSettingsV220,
    initial_gaussian_spinor_v220,
    run_spinor_exact_grid_v220,
)
from .v220_benchmark import run_v0220_release_benchmark


@dataclass(frozen=True)
class V221AcceptanceThresholds:
    maximum_component_derivative_error: float = 2.0e-9
    maximum_fine_basis_population_difference: float = 1.0e-8
    maximum_basis_narrowing_ratio: float = 0.05
    maximum_fine_sparse_coefficient_error: float = 1.0e-12
    monotonic_tolerance: float = 1.0e-14


class _GenericThreeStateTwoCoordinateProviderV221:
    """No-config provider proving that the audit derives dimensions from data."""

    def __init__(self):
        self._projectors = {
            f"root_{index}": np.diag(
                [1.0 if row == index else 0.0 for row in range(3)]
            ).astype(complex)
            for index in range(3)
        }
        self._symmetry = SOCSymmetryContractV221(
            "even", np.eye(3, dtype=complex), self._projectors
        )
        model_space = ElectronicModelSpaceV213(
            name="v0.22.1 generic three-singlet two-coordinate probe",
            representation="fixed_general",
            states=tuple(
                ElectronicStateDescriptorV213(
                    f"S{index}", f"S{index}", 1, "M=0", 0
                )
                for index in range(3)
            ),
            complete_multiplets=True,
        ).validate()
        self.provenance = ElectronicOperatorProvenanceV213(
            model_name="v0.22.1 generic dimension-neutral SOC probe",
            model_version="1",
            model_space=model_space,
            spin_free_method="analytic real linear diagonal surfaces",
            soc_enabled=True,
            soc_method="analytic real linear symmetric SOC",
            scalar_relativistic_method="none",
            derivative_method="analytic component derivatives",
            parameters={
                **self._symmetry.as_provenance_parameters(),
                "soc_signal_expected": True,
                "basis_order": ["S0", "S1", "S2"],
            },
        ).validate()

    @property
    def time_reversal_matrix(self):
        return self._symmetry.time_reversal_matrix.copy()

    @property
    def projectors(self):
        return {name: value.copy() for name, value in self._projectors.items()}

    @property
    def soc_symmetry_contract(self):
        return self._symmetry

    def components(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (2,) or not np.all(np.isfinite(q)):
            raise ValueError("generic SOC probe requires two finite coordinates.")
        q0, q1 = q
        H0 = np.diag(
            [
                0.01 + 0.012 * q0 - 0.004 * q1,
                0.02 - 0.007 * q0 + 0.006 * q1,
                0.03 + 0.005 * q0 + 0.009 * q1,
            ]
        ).astype(complex)
        K0 = np.asarray(
            [
                np.diag([0.012, -0.007, 0.005]),
                np.diag([-0.004, 0.006, 0.009]),
            ],
            dtype=complex,
        )
        Hso = np.zeros((3, 3), dtype=complex)
        Kso = np.zeros((2, 3, 3), dtype=complex)
        values = (
            (0, 1, 0.0010 + 0.0002 * q0 - 0.0001 * q1, 0.0002, -0.0001),
            (0, 2, -0.0006 + 0.0001 * q0 + 0.0003 * q1, 0.0001, 0.0003),
            (1, 2, 0.0008 - 0.0002 * q0 + 0.0002 * q1, -0.0002, 0.0002),
        )
        for left, right, value, derivative_0, derivative_1 in values:
            Hso[left, right] = Hso[right, left] = value
            Kso[0, left, right] = Kso[0, right, left] = derivative_0
            Kso[1, left, right] = Kso[1, right, left] = derivative_1
        return SOCOperatorComponentsV220(q, H0, K0, Hso, Kso).validate()

    def evaluate_snapshot(self, q):
        components = self.components(q)
        point = compose_electronic_operator_v213(
            q=components.q,
            H_spin_free=components.H_spin_free,
            dH_spin_free_dq=components.K_spin_free,
            H_soc=components.H_soc,
            dH_soc_dq=components.K_soc,
            connection_q=np.zeros((2, 3, 3), dtype=complex),
            mass_matrix_q_au=np.diag([800.0, 900.0]),
            provenance=self.provenance,
        )
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=np.eye(3, dtype=complex),
            metadata={"provider": "_GenericThreeStateTwoCoordinateProviderV221"},
        ).validate()

    def snapshot_overlap(self, left, right):
        return np.eye(3, dtype=complex)


class _CancelledComponentDerivativeV221:
    def __init__(self, base):
        self.base = base
        self.provenance = base.provenance

    @property
    def time_reversal_matrix(self):
        return self.base.time_reversal_matrix

    @property
    def projectors(self):
        return self.base.projectors

    @property
    def soc_symmetry_contract(self):
        return self.base.soc_symmetry_contract

    def components(self, q):
        components = self.base.components(q)
        indices = np.arange(1, 5, dtype=float)
        vector = (1.0 + 0.13 * indices) + 1j * (0.19 - 0.07 * indices)
        vector /= np.linalg.norm(vector)
        weights = abs(vector) ** 2
        delta = 1.0e-3 * np.diag(
            [0.0, weights[2], -(weights[1] + weights[3]), weights[2]]
        )
        if abs(np.vdot(vector, delta @ vector)) > 1.0e-16:
            raise AssertionError("cancelled derivative fixture is not force-orthogonal.")
        return SOCOperatorComponentsV220(
            components.q,
            components.H_spin_free,
            components.K_spin_free - delta[None, :, :],
            components.H_soc,
            components.K_soc + delta[None, :, :],
        ).validate()

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left, right)


class _NoConfigProviderV221:
    def __init__(self, base):
        self.base = base
        self.provenance = base.provenance

    @property
    def time_reversal_matrix(self):
        return self.base.time_reversal_matrix

    @property
    def projectors(self):
        return self.base.projectors

    @property
    def soc_symmetry_contract(self):
        return self.base.soc_symmetry_contract

    def components(self, q):
        return self.base.components(q)

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left, right)


def _operator_hardening_v221():
    generic = _GenericThreeStateTwoCoordinateProviderV221()
    generic_report = audit_physical_soc_provider_v220(
        generic, np.asarray([0.13, -0.21]), fermionic=False
    )
    cancelled = audit_physical_soc_provider_v220(
        _CancelledComponentDerivativeV221(
            AnalyticSingletTripletSOCProviderV220()
        ),
        np.asarray([0.17]),
        fermionic=False,
    )
    no_config = audit_physical_soc_provider_v220(
        _NoConfigProviderV221(AnalyticDoubletSOCProviderV220()),
        np.asarray([-0.11]),
        fermionic=True,
    )
    return {
        "generic_three_state_two_coordinate": generic_report.as_dict(),
        "cancelled_component_derivative": cancelled.as_dict(),
        "no_config_provider": no_config.as_dict(),
    }


def _symmetry_hardening_v221():
    mixed_space = ElectronicModelSpaceV213(
        name="mixed electron parity negative fixture",
        representation="fixed_spin_diabatic",
        states=(
            ElectronicStateDescriptorV213("S", "S", 1, "M=0", 0),
            ElectronicStateDescriptorV213("D+", "D", 2, "M=+1/2", 0),
            ElectronicStateDescriptorV213("D-", "D", 2, "M=-1/2", 0),
        ),
        complete_multiplets=True,
    ).validate()
    mixed_contract = SOCSymmetryContractV221(
        "even", np.eye(3, dtype=complex), {"all": np.eye(3, dtype=complex)}
    )
    mixed = audit_soc_symmetry_contract_v221(mixed_space, mixed_contract)

    two_singlets = ElectronicModelSpaceV213(
        name="nonunitary time-reversal negative fixture",
        representation="fixed_general",
        states=(
            ElectronicStateDescriptorV213("S0", "S0", 1, "M=0", 0),
            ElectronicStateDescriptorV213("S1", "S1", 1, "M=0", 0),
        ),
        complete_multiplets=True,
    ).validate()
    bad_J = np.asarray([[1.0, 2.0], [0.0, -1.0]], dtype=complex)
    nonunitary = audit_soc_symmetry_contract_v221(
        two_singlets,
        SOCSymmetryContractV221(
            "even", bad_J, {"all": np.eye(2, dtype=complex)}
        ),
    )

    provider = AnalyticSingletTripletSOCProviderV220()
    parameters = copy.deepcopy(provider.provenance.parameters)
    parameters["physical_projectors"] = {
        "singlet": provider.projectors["triplet"].tolist(),
        "triplet": provider.projectors["singlet"].tolist(),
    }
    mismatched_provenance = replace(provider.provenance, parameters=parameters)
    mismatch = audit_soc_symmetry_contract_v221(
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
        provenance=mismatched_provenance,
        fermionic=False,
    )
    valid = audit_soc_symmetry_contract_v221(
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
        provenance=provider.provenance,
        fermionic=False,
    )
    return {
        "mixed_parity": mixed.as_dict(),
        "nonunitary_time_reversal": nonunitary.as_dict(),
        "provenance_mismatch": mismatch.as_dict(),
        "valid_reference": valid.as_dict(),
    }


def _grid_hardening_v221():
    x = np.linspace(-4.0, 4.0, 64, endpoint=False)
    psi0 = initial_gaussian_spinor_v220(
        x, np.asarray([1.0, 0.0, 0.0, 0.0])
    )
    no_config_provider = _NoConfigProviderV221(
        AnalyticSingletTripletSOCProviderV220()
    )
    output = run_spinor_exact_grid_v220(
        no_config_provider,
        x,
        psi0,
        settings=SpinorGridSettingsV220(dt=0.01, steps=5, store_every=2),
    )

    base = AnalyticDoubletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, 22107),
        np.asarray([[0.11], [-0.08], [0.17], [-0.13]]),
        np.asarray([0.20, -0.31, 0.14, -0.09]),
    )
    local_provenance = DoubletSOCConfigV220().provenance("local_general")
    moving = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(base, gauge), local_provenance
    )
    moving_rejected = False
    try:
        run_spinor_exact_grid_v220(
            moving,
            x,
            psi0,
            settings=SpinorGridSettingsV220(dt=0.01, steps=1, store_every=1),
        )
    except ValueError as exc:
        moving_rejected = "fixed electronic frame" in str(exc)
    return {
        "requested_final_time": 0.05,
        "recorded_final_time": float(output["time"][-1]),
        "recorded_times": output["time"].tolist(),
        "provider_has_config": hasattr(no_config_provider, "config"),
        "provider_mass_fingerprint": output["provider_fingerprint"],
        "fixed_frame_certified": bool(output["fixed_frame_certified"]),
        "constant_mass_certified": bool(output["constant_mass_certified"]),
        "moving_frame_rejected": moving_rejected,
    }


def _block_projector_population_v221(coefficients, metric, projector, nbasis):
    coefficients = np.asarray(coefficients, dtype=complex)
    metric = metric.toarray() if hasattr(metric, "toarray") else np.asarray(metric)
    block_projector = np.kron(np.eye(int(nbasis)), np.asarray(projector, complex))
    numerator = np.real(
        np.vdot(coefficients, metric @ (block_projector @ coefficients))
    )
    denominator = np.real(np.vdot(coefficients, metric @ coefficients))
    return float(numerator / denominator)


def _initial_projection_v221(provider, basis, *, nx=512):
    x = np.linspace(-8.0, 8.0, int(nx), endpoint=False)
    dx = float(x[1] - x[0])
    psi = initial_gaussian_spinor_v220(
        x,
        np.asarray([1.0, 0.0, 0.0, 0.0]),
        center=-1.0,
        momentum=1.2,
        width=0.7,
    )
    projection = project_grid_wavefunction_fixed_frame_v213(
        psi.T, x[:, None], dx, basis
    )
    return projection


def _basis_convergence_v221():
    provider = AnalyticSingletTripletSOCProviderV220()
    ladders = (
        (0.0,),
        (-0.7, 0.0, 0.7),
        (-1.4, -0.7, 0.0, 0.7, 1.4),
    )
    rows = []
    for ladder_index, shifts in enumerate(ladders):
        basis = [
            BlockMolecularTBFV21(
                1000 + 10 * ladder_index + index,
                np.asarray([-1.0 + shift]),
                np.asarray([1.2]),
                np.asarray([[0.7]]),
            )
            for index, shift in enumerate(shifts)
        ]
        projection = _initial_projection_v221(provider, basis)
        output = run_prescribed_block_dynamics_v21(
            basis,
            projection.coefficients,
            provider,
            np.zeros((len(basis), 1)),
            dt=0.01,
            steps=100,
            store_every=100,
            settings=PrescribedBlockDynamicsSettingsV21(use_dense_reference=True),
        )
        population = _block_projector_population_v221(
            output["final_coefficients"],
            output["final_S"],
            provider.projectors["triplet"],
            len(basis),
        )
        rows.append(
            {
                "basis_size": len(basis),
                "initial_projection_fidelity": projection.fidelity,
                "final_triplet_population": population,
            }
        )
    coarse_difference = abs(
        rows[0]["final_triplet_population"]
        - rows[1]["final_triplet_population"]
    )
    fine_difference = abs(
        rows[1]["final_triplet_population"]
        - rows[2]["final_triplet_population"]
    )
    return {
        "rows": rows,
        "coarse_population_difference": coarse_difference,
        "fine_population_difference": fine_difference,
        "narrowing_ratio": fine_difference / max(coarse_difference, 1.0e-30),
    }


def _metric_coefficient_error_v221(reference, candidate, metric):
    metric = metric.toarray() if hasattr(metric, "toarray") else np.asarray(metric)
    overlap = np.vdot(reference, metric @ candidate)
    phase = 1.0 if abs(overlap) < 1.0e-30 else np.exp(-1j * np.angle(overlap))
    difference = phase * candidate - reference
    numerator = max(float(np.real(np.vdot(difference, metric @ difference))), 0.0)
    denominator = max(float(np.real(np.vdot(reference, metric @ reference))), 1.0e-30)
    return float(np.sqrt(numerator / denominator))


def _sparse_threshold_convergence_v221():
    provider = AnalyticDoubletSOCProviderV220()
    shifts = (-1.6, -0.8, 0.0, 0.8, 1.6, 2.4)
    basis = [
        BlockMolecularTBFV21(
            2000 + index,
            np.asarray([-1.0 + shift]),
            np.asarray([1.2]),
            np.asarray([[0.7]]),
        )
        for index, shift in enumerate(shifts)
    ]
    projection = _initial_projection_v221(provider, basis)
    velocities = np.zeros((len(basis), 1))
    common = dict(dt=0.01, steps=20, store_every=20)
    dense = run_prescribed_block_dynamics_v21(
        basis,
        projection.coefficients,
        provider,
        velocities,
        settings=PrescribedBlockDynamicsSettingsV21(use_dense_reference=True),
        **common,
    )
    rows = []
    for threshold in (1.2, 0.5, 0.15, 0.05):
        graph = BlockSparseSettingsV21(
            enter_score=threshold,
            exit_score=0.7 * threshold,
            search_overlap_floor=1.0e-10,
            local_omitted_score_l2_budget=1.0e9,
            use_kdtree=False,
        )
        output = run_prescribed_block_dynamics_v21(
            basis,
            projection.coefficients,
            provider,
            velocities,
            settings=PrescribedBlockDynamicsSettingsV21(
                graph=graph, use_dense_reference=False
            ),
            **common,
        )
        rows.append(
            {
                "enter_score": threshold,
                "active_edges": len(output["final_active_edges"]),
                "coefficient_error": _metric_coefficient_error_v221(
                    dense["final_coefficients"],
                    output["final_coefficients"],
                    dense["final_S"],
                ),
            }
        )
    return {"rows": rows}


def run_v0221_release_benchmark():
    inherited = run_v0220_release_benchmark()
    operator = _operator_hardening_v221()
    symmetry = _symmetry_hardening_v221()
    grid = _grid_hardening_v221()
    basis = _basis_convergence_v221()
    sparse = _sparse_threshold_convergence_v221()
    thresholds = V221AcceptanceThresholds()
    generic = operator["generic_three_state_two_coordinate"]
    cancelled = operator["cancelled_component_derivative"]
    no_config = operator["no_config_provider"]
    sparse_errors = [row["coefficient_error"] for row in sparse["rows"]]
    sparse_edges = [row["active_edges"] for row in sparse["rows"]]

    new_checks = {
        "derivatives::full_matrix_spin_free_components": generic[
            "maximum_spin_free_component_derivative_error"
        ] <= thresholds.maximum_component_derivative_error,
        "derivatives::full_matrix_SOC_components": generic[
            "maximum_soc_component_derivative_error"
        ] <= thresholds.maximum_component_derivative_error,
        "derivatives::cancelled_component_error_rejected": (
            not cancelled["passed"]
            and cancelled["checks"]["H_decomposition"]
            and cancelled["checks"]["K_decomposition"]
            and cancelled["checks"]["cross_geometry_differentials"]
            and cancelled["checks"]["SOC_force_derivative"]
            and not cancelled["checks"]["spin_free_component_derivatives"]
            and not cancelled["checks"]["SOC_component_derivatives"]
        ),
        "contract::arbitrary_state_and_coordinate_dimensions": (
            generic["passed"] and len(generic["q"]) == 2
        ),
        "contract::provider_config_independence": no_config["passed"],
        "symmetry::mixed_electron_parity_rejected": (
            not symmetry["mixed_parity"]["passed"]
            and not symmetry["mixed_parity"]["checks"]["single_electron_parity"]
        ),
        "symmetry::time_reversal_unitarity_required": (
            not symmetry["nonunitary_time_reversal"]["passed"]
            and symmetry["nonunitary_time_reversal"][
                "time_reversal_square_residual"
            ] == 0.0
            and not symmetry["nonunitary_time_reversal"]["checks"][
                "time_reversal_unitarity"
            ]
        ),
        "symmetry::numerical_provenance_identity": (
            symmetry["valid_reference"]["passed"]
            and not symmetry["provenance_mismatch"]["passed"]
            and not symmetry["provenance_mismatch"]["checks"][
                "symmetry_provenance_identity"
            ]
        ),
        "grid::final_step_always_recorded": (
            grid["recorded_final_time"] == grid["requested_final_time"]
            and grid["recorded_times"] == [0.0, 0.02, 0.04, 0.05]
        ),
        "grid::moving_frame_rejected": grid["moving_frame_rejected"],
        "grid::contract_mass_without_config": (
            not grid["provider_has_config"]
            and grid["fixed_frame_certified"]
            and grid["constant_mass_certified"]
        ),
        "convergence::SOC_Gaussian_basis": (
            basis["fine_population_difference"]
            <= thresholds.maximum_fine_basis_population_difference
            and basis["narrowing_ratio"]
            <= thresholds.maximum_basis_narrowing_ratio
        ),
        "convergence::SOC_sparse_threshold": (
            all(
                later <= earlier + thresholds.monotonic_tolerance
                for earlier, later in zip(sparse_errors[:-1], sparse_errors[1:])
            )
            and all(
                later >= earlier
                for earlier, later in zip(sparse_edges[:-1], sparse_edges[1:])
            )
            and sparse_errors[-1]
            <= thresholds.maximum_fine_sparse_coefficient_error
        ),
        "campaign::native_boolean_checks": all(
            type(value) is bool
            for value in inherited["acceptance"]["checks"].values()
        ),
    }
    new_checks = {name: bool(value) for name, value in new_checks.items()}
    if len(new_checks) != 14:
        raise AssertionError("v0.22.1 campaign must define exactly 14 new gates.")
    inherited_checks = {
        f"inherited_v0220::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    checks = {**inherited_checks, **new_checks}
    if len(checks) != 67:
        raise AssertionError("v0.22.1 campaign must define exactly 67 total gates.")
    return {
        "release": "v0.22.1",
        "theme": "corrective SOC derivative, symmetry, grid, and convergence hardening",
        "operator_hardening": operator,
        "symmetry_hardening": symmetry,
        "grid_hardening": grid,
        "gaussian_basis_convergence": basis,
        "sparse_threshold_convergence": sparse,
        "inherited_v0220": inherited,
        "soc": {
            "physical_analytic_SOC": True,
            "analytic_models_only": True,
            "ab_initio_SOC_validated": False,
            "molecular_SOC_backend_admitted": False,
            "spin_free_mode_permanent": True,
        },
        "pyscf": {
            "installed_in_build_environment": bool(importlib.util.find_spec("pyscf")),
            "runtime_validated": False,
        },
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
