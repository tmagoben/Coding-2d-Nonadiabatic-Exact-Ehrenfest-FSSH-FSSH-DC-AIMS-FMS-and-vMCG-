"""Canonical 53-gate physical analytic-SOC campaign for v0.22.0."""

from dataclasses import asdict, dataclass
import importlib.util
from pathlib import Path
import tempfile
import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
    SOCOperatorComponentsV220,
    SingletTripletSOCConfigV220,
)
from .block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,
    BlockSparseSettingsV21,
)
from .checkpoint_restart_v214 import (
    SelfConsistentBlockSettingsV214,
    load_self_consistent_checkpoint_v214,
    run_self_consistent_block_dynamics_v214,
    save_self_consistent_checkpoint_v214,
)
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from .electronic_contract_v213 import ContractedElectronicOperatorProviderV213
from .physical_soc_validation_v220 import (
    audit_kramers_degeneracy_v220,
    audit_physical_soc_provider_v220,
    projector_population_v220,
    time_reversal_residual_v220,
    transform_projector_v220,
    transform_time_reversal_matrix_v220,
)
from .provider_differential_audit_v214 import audit_provider_differentials_v214
from .spinor_exact_grid_v220 import (
    SpinorGridSettingsV220,
    initial_gaussian_spinor_v220,
    phase_aligned_spinor_grid_error_v220,
    run_spinor_exact_grid_v220,
)
from .v214_benchmark import run_v0214_release_benchmark


@dataclass(frozen=True)
class V220AcceptanceThresholds:
    maximum_composition_error: float = 1.0e-13
    maximum_time_reversal_residual: float = 1.0e-12
    maximum_force_error: float = 2.0e-10
    maximum_kramers_splitting: float = 1.0e-11
    maximum_grid_norm_drift: float = 5.0e-14
    maximum_grid_energy_drift: float = 1.0e-12
    minimum_timestep_order: float = 1.99
    maximum_grid_population_resolution_error: float = 1.0e-12
    maximum_gaussian_grid_population_error: float = 1.0e-8
    maximum_restart_error: float = 2.0e-12
    maximum_gauge_error: float = 2.0e-12


def _basis_v220():
    return [
        BlockMolecularTBFV21(
            3, np.asarray([-0.55]), np.asarray([0.20]), np.asarray([[1.15]])
        ),
        BlockMolecularTBFV21(
            8, np.asarray([0.45]), np.asarray([-0.10]), np.asarray([[1.35]])
        ),
    ]


def _coefficients_v220():
    return np.asarray(
        [0.76 + 0.08j, 0.0, 0.0, 0.0, 0.22 - 0.11j, 0.0, 0.0, 0.0]
    )


def _dense_settings_v220():
    return SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )


def _phase_aligned_metric_error_v220(reference, candidate, metric):
    overlap = np.vdot(reference, metric @ candidate)
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else np.exp(-1j * np.angle(overlap))
    difference = phase * candidate - reference
    return float(np.sqrt(max(np.real(np.vdot(difference, metric @ difference)), 0.0)))


def _trajectory_errors_v220(reference, candidate):
    return {
        "position": float(
            max(
                np.linalg.norm(left.q - right.q)
                for left, right in zip(reference["final_basis"], candidate["final_basis"])
            )
        ),
        "momentum": float(
            max(
                np.linalg.norm(left.p - right.p)
                for left, right in zip(reference["final_basis"], candidate["final_basis"])
            )
        ),
        "coefficient": _phase_aligned_metric_error_v220(
            reference["final_coefficients"],
            candidate["final_coefficients"],
            reference["final_S"],
        ),
    }


class _WrongSOCDerivativeV220:
    def __init__(self, base):
        self.base = base
        self.config = base.config
        self.provenance = base.provenance

    @property
    def time_reversal_matrix(self):
        return self.base.time_reversal_matrix

    @property
    def projectors(self):
        return self.base.projectors

    def components(self, q):
        components = self.base.components(q)
        wrong = components.K_soc.copy()
        wrong[0] += 1.0e-3 * np.eye(4)
        return SOCOperatorComponentsV220(
            components.q,
            components.H_spin_free,
            components.K_spin_free,
            components.H_soc,
            wrong,
        ).validate()

    def evaluate_snapshot(self, q):
        return self.base.evaluate_snapshot(q)

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left, right)


def _gauge_v220():
    return PhaseMixingGaugeV21(
        random_unitary_v21(4, 22031),
        np.asarray([[0.11], [-0.08], [0.17], [-0.13]]),
        np.asarray([0.20, -0.31, 0.14, -0.09]),
    )


def _operator_campaign_v220():
    singlet_triplet = AnalyticSingletTripletSOCProviderV220()
    doublet = AnalyticDoubletSOCProviderV220()
    st_audit = audit_physical_soc_provider_v220(
        singlet_triplet, np.asarray([0.17]), fermionic=False
    )
    doublet_audit = audit_physical_soc_provider_v220(
        doublet, np.asarray([-0.23]), fermionic=True
    )
    kramers = audit_kramers_degeneracy_v220(
        doublet,
        [np.asarray([-1.2]), np.asarray([0.0]), np.asarray([1.1])],
    )
    constant_provider = AnalyticSingletTripletSOCProviderV220(
        SingletTripletSOCConfigV220(
            lambda_real_gradient=0.0,
            lambda_imag_gradient=0.0,
            lambda_zero_gradient=0.0,
        )
    )
    constant_components = constant_provider.components(np.asarray([0.31]))
    zero_equivalence = {}
    for name, config_type, provider_type in (
        (
            "singlet_triplet",
            SingletTripletSOCConfigV220,
            AnalyticSingletTripletSOCProviderV220,
        ),
        ("doublet", DoubletSOCConfigV220, AnalyticDoubletSOCProviderV220),
    ):
        enabled = provider_type(config_type(soc_scale=0.0, soc_enabled=True))
        disabled = provider_type(config_type(soc_scale=0.0, soc_enabled=False))
        q = np.asarray([0.29])
        left = enabled.evaluate_snapshot(q).point
        right = disabled.evaluate_snapshot(q).point
        zero_equivalence[name] = {
            "H_error": float(np.max(np.abs(left.H - right.H))),
            "K_error": float(np.max(np.abs(left.dH_dq - right.dH_dq))),
            "D_error": float(np.max(np.abs(left.connection_q - right.connection_q))),
            "mass_error": float(
                np.max(np.abs(left.mass_matrix_q_au - right.mass_matrix_q_au))
            ),
        }
    wrong = audit_physical_soc_provider_v220(
        _WrongSOCDerivativeV220(AnalyticSingletTripletSOCProviderV220()),
        np.asarray([0.17]),
        fermionic=False,
    )
    broken_H = doublet.evaluate_snapshot(np.asarray([0.1])).point.H.copy()
    broken_H[0, 0] += 1.0e-3
    broken_residual = time_reversal_residual_v220(
        broken_H, doublet.time_reversal_matrix
    )
    broken_energies = np.linalg.eigvalsh(broken_H)
    broken_splitting = float(
        max(abs(broken_energies[1::2] - broken_energies[0::2]))
    )
    return {
        "singlet_triplet": st_audit.as_dict(),
        "doublet": doublet_audit.as_dict(),
        "kramers": kramers.as_dict(),
        "constant_SOC": {
            "H_soc_norm": float(np.linalg.norm(constant_components.H_soc)),
            "K_soc_norm": float(np.linalg.norm(constant_components.K_soc)),
        },
        "zero_SOC_equivalence": zero_equivalence,
        "wrong_SOC_derivative": wrong.as_dict(),
        "broken_kramers": {
            "time_reversal_residual": broken_residual,
            "pair_splitting": broken_splitting,
        },
    }


def _gauge_campaign_v220():
    config = DoubletSOCConfigV220()
    base = AnalyticDoubletSOCProviderV220(config)
    gauge = _gauge_v220()
    q = np.asarray([0.19])
    G = gauge.matrix(q)
    point = base.evaluate_snapshot(q).point
    transformed_H = G.conj().T @ point.H @ G
    transformed_J = transform_time_reversal_matrix_v220(
        base.time_reversal_matrix, G
    )
    transformed_J_residual = time_reversal_residual_v220(
        transformed_H, transformed_J
    )
    untransformed_J_residual = time_reversal_residual_v220(
        transformed_H, base.time_reversal_matrix
    )
    vector = np.asarray([0.55 + 0.1j, -0.12j, 0.31 - 0.2j, 0.42 + 0.08j])
    transformed_vector = G.conj().T @ vector
    projector_errors = []
    for projector in base.projectors.values():
        projector_errors.append(
            abs(
                projector_population_v220(vector, projector)
                - projector_population_v220(
                    transformed_vector, transform_projector_v220(projector, G)
                )
            )
        )
    local_provenance = config.provenance("local_general")
    local = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(
            AnalyticDoubletSOCProviderV220(config), gauge
        ),
        local_provenance,
    )
    differential = audit_provider_differentials_v214(
        local, q, local_provenance
    )

    basis = _basis_v220()
    coefficients = _coefficients_v220()
    local_coefficients = np.concatenate(
        [
            gauge.matrix(item.q).conj().T
            @ coefficients[4 * index : 4 * index + 4]
            for index, item in enumerate(basis)
        ]
    )
    common = dict(
        dt=0.002,
        steps=6,
        store_every=2,
        settings=_dense_settings_v220(),
    )
    fixed_run = run_self_consistent_block_dynamics_v214(
        base,
        base.provenance,
        initial_basis=basis,
        C0=coefficients,
        **common,
    )
    local_run = run_self_consistent_block_dynamics_v214(
        local,
        local_provenance,
        initial_basis=basis,
        C0=local_coefficients,
        **common,
    )
    mapped = np.concatenate(
        [
            gauge.matrix(item.q)
            @ local_run["final_coefficients"][4 * index : 4 * index + 4]
            for index, item in enumerate(local_run["final_basis"])
        ]
    )
    dynamics_errors = {
        "position": float(
            max(
                np.linalg.norm(left.q - right.q)
                for left, right in zip(fixed_run["final_basis"], local_run["final_basis"])
            )
        ),
        "momentum": float(
            max(
                np.linalg.norm(left.p - right.p)
                for left, right in zip(fixed_run["final_basis"], local_run["final_basis"])
            )
        ),
        "coefficient": _phase_aligned_metric_error_v220(
            fixed_run["final_coefficients"], mapped, fixed_run["final_S"]
        ),
    }
    return {
        "transformed_time_reversal_residual": transformed_J_residual,
        "untransformed_time_reversal_residual": untransformed_J_residual,
        "maximum_projector_population_error": max(projector_errors, default=0.0),
        "differential": differential.as_dict(),
        "dynamics_errors": dynamics_errors,
    }


def _short_gaussian_population_v220(provider_type, projector_name):
    provider = provider_type()
    basis = [
        BlockMolecularTBFV21(
            1, np.asarray([-1.0]), np.asarray([1.2]), np.asarray([[0.7]])
        )
    ]
    output = run_self_consistent_block_dynamics_v214(
        provider,
        provider.provenance,
        initial_basis=basis,
        C0=np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex),
        dt=0.002,
        steps=100,
        store_every=100,
        settings=SelfConsistentBlockSettingsV214(
            use_dense_reference=True,
            corrector_iterations=4,
            momentum_tolerance=1.0e-12,
        ),
    )
    coefficients = output["final_coefficients"]
    projector = provider.projectors[projector_name]
    return float(
        np.real(np.vdot(coefficients, projector @ coefficients))
        / np.real(np.vdot(coefficients, coefficients))
    )


def _grid_run_v220(provider, x, *, dt=0.04, final_time=4.0):
    psi0 = initial_gaussian_spinor_v220(
        x,
        np.asarray([1.0, 0.0, 0.0, 0.0]),
        center=-1.0,
        momentum=1.2,
        width=0.7,
    )
    steps = int(round(float(final_time) / float(dt)))
    return run_spinor_exact_grid_v220(
        provider,
        x,
        psi0,
        settings=SpinorGridSettingsV220(dt=dt, steps=steps, store_every=steps),
    )


def _exact_grid_campaign_v220():
    x = np.linspace(-8.0, 8.0, 256, endpoint=False)
    st = _grid_run_v220(AnalyticSingletTripletSOCProviderV220(), x)
    doublet = _grid_run_v220(AnalyticDoubletSOCProviderV220(), x)
    timestep_outputs = [
        _grid_run_v220(AnalyticDoubletSOCProviderV220(), x, dt=dt)
        for dt in (0.08, 0.04, 0.02)
    ]
    dx = x[1] - x[0]
    coarse_difference = phase_aligned_spinor_grid_error_v220(
        timestep_outputs[0]["psi"][-1], timestep_outputs[1]["psi"][-1], dx
    )
    fine_difference = phase_aligned_spinor_grid_error_v220(
        timestep_outputs[1]["psi"][-1], timestep_outputs[2]["psi"][-1], dx
    )
    order = float(np.log(coarse_difference / fine_difference) / np.log(2.0))
    fine_x = np.linspace(-8.0, 8.0, 512, endpoint=False)
    wide_x = np.linspace(-10.0, 10.0, 320, endpoint=False)
    fine = _grid_run_v220(AnalyticDoubletSOCProviderV220(), fine_x)
    wide = _grid_run_v220(AnalyticDoubletSOCProviderV220(), wide_x)
    base_population = float(doublet["populations"]["doublet_2"][-1])
    resolution_error = max(
        abs(base_population - float(fine["populations"]["doublet_2"][-1])),
        abs(base_population - float(wide["populations"]["doublet_2"][-1])),
    )
    short_time = 0.2
    short_settings = SpinorGridSettingsV220(dt=0.01, steps=20, store_every=20)
    gaussian_grid_errors = {}
    for name, provider_type, projector_name in (
        (
            "singlet_triplet",
            AnalyticSingletTripletSOCProviderV220,
            "triplet",
        ),
        ("doublet", AnalyticDoubletSOCProviderV220, "doublet_2"),
    ):
        exact_x = np.linspace(-8.0, 8.0, 512, endpoint=False)
        exact_psi = initial_gaussian_spinor_v220(
            exact_x,
            np.asarray([1.0, 0.0, 0.0, 0.0]),
            center=-1.0,
            momentum=1.2,
            width=0.7,
        )
        exact = run_spinor_exact_grid_v220(
            provider_type(), exact_x, exact_psi, settings=short_settings
        )
        exact_population = float(exact["populations"][projector_name][-1])
        gaussian_population = _short_gaussian_population_v220(
            provider_type, projector_name
        )
        gaussian_grid_errors[name] = {
            "time": short_time,
            "gaussian_population": gaussian_population,
            "exact_grid_population": exact_population,
            "absolute_error": abs(gaussian_population - exact_population),
        }
    return {
        "singlet_triplet": {
            "maximum_norm_drift": st["maximum_norm_drift"],
            "maximum_energy_drift": st["maximum_energy_drift"],
            "final_populations": {
                name: float(values[-1]) for name, values in st["populations"].items()
            },
        },
        "doublet": {
            "maximum_norm_drift": doublet["maximum_norm_drift"],
            "maximum_energy_drift": doublet["maximum_energy_drift"],
            "final_populations": {
                name: float(values[-1])
                for name, values in doublet["populations"].items()
            },
        },
        "timestep_convergence": {
            "coarse_difference": coarse_difference,
            "fine_difference": fine_difference,
            "observed_order": order,
        },
        "grid_resolution_and_box_population_error": resolution_error,
        "gaussian_exact_grid_population_errors": gaussian_grid_errors,
    }


def _restart_campaign_v220():
    provider = AnalyticSingletTripletSOCProviderV220()
    common = dict(dt=0.002, store_every=2, settings=_dense_settings_v220())
    full = run_self_consistent_block_dynamics_v214(
        provider,
        provider.provenance,
        initial_basis=_basis_v220(),
        C0=_coefficients_v220(),
        steps=8,
        **common,
    )
    first_provider = AnalyticSingletTripletSOCProviderV220()
    first = run_self_consistent_block_dynamics_v214(
        first_provider,
        first_provider.provenance,
        initial_basis=_basis_v220(),
        C0=_coefficients_v220(),
        steps=3,
        **common,
    )
    resumed_provider = AnalyticSingletTripletSOCProviderV220()
    resumed = run_self_consistent_block_dynamics_v214(
        resumed_provider,
        resumed_provider.provenance,
        checkpoint=first["checkpoint"],
        steps=5,
        store_every=1,
        settings=_dense_settings_v220(),
    )

    sparse_control = BlockSparseSettingsV21(
        enter_score=0.30,
        exit_score=0.10,
        search_overlap_floor=1.0e-8,
        local_omitted_score_l2_budget=1.0,
        use_kdtree=False,
    )
    sparse_settings = SelfConsistentBlockSettingsV214(
        graph=sparse_control,
        use_dense_reference=False,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    sparse_provider = AnalyticDoubletSOCProviderV220()
    sparse_full = run_self_consistent_block_dynamics_v214(
        sparse_provider,
        sparse_provider.provenance,
        initial_basis=_basis_v220(),
        C0=_coefficients_v220(),
        dt=0.002,
        steps=6,
        store_every=2,
        settings=sparse_settings,
    )
    sparse_first_provider = AnalyticDoubletSOCProviderV220()
    sparse_first = run_self_consistent_block_dynamics_v214(
        sparse_first_provider,
        sparse_first_provider.provenance,
        initial_basis=_basis_v220(),
        C0=_coefficients_v220(),
        dt=0.002,
        steps=2,
        store_every=1,
        settings=sparse_settings,
    )
    sparse_resumed_provider = AnalyticDoubletSOCProviderV220()
    sparse_resumed = run_self_consistent_block_dynamics_v214(
        sparse_resumed_provider,
        sparse_resumed_provider.provenance,
        checkpoint=sparse_first["checkpoint"],
        steps=4,
        store_every=2,
        settings=sparse_settings,
    )

    changed_rejected = False
    changed = AnalyticSingletTripletSOCProviderV220(
        SingletTripletSOCConfigV220(soc_scale=0.0031)
    )
    try:
        run_self_consistent_block_dynamics_v214(
            changed,
            changed.provenance,
            checkpoint=first["checkpoint"],
            steps=1,
            settings=_dense_settings_v220(),
        )
    except ValueError:
        changed_rejected = True

    with tempfile.TemporaryDirectory(prefix="v220-soc-checkpoint-") as directory:
        path = save_self_consistent_checkpoint_v214(
            Path(directory) / "soc.npz", first["checkpoint"]
        )
        with np.load(path, allow_pickle=False) as archive:
            arrays = {name: archive[name].copy() for name in archive.files}
        arrays["q"][0, 0] += 1.0e-3
        corrupted = Path(directory) / "soc-corrupted.npz"
        with corrupted.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
        corruption_rejected = False
        try:
            load_self_consistent_checkpoint_v214(corrupted)
        except ValueError:
            corruption_rejected = True
    return {
        "dense_errors": _trajectory_errors_v220(full, resumed),
        "sparse_errors": _trajectory_errors_v220(sparse_full, sparse_resumed),
        "sparse_edges": {
            "checkpoint": sparse_first["checkpoint"].active_uid_edges.tolist(),
            "uninterrupted": [list(edge) for edge in sparse_full["final_active_uid_edges"]],
            "resumed": [list(edge) for edge in sparse_resumed["final_active_uid_edges"]],
        },
        "changed_soc_provenance_rejected": changed_rejected,
        "corruption_rejected": corruption_rejected,
    }


def run_v0220_release_benchmark():
    inherited = run_v0214_release_benchmark()
    operator = _operator_campaign_v220()
    gauge = _gauge_campaign_v220()
    exact = _exact_grid_campaign_v220()
    restart = _restart_campaign_v220()
    thresholds = V220AcceptanceThresholds()
    st = operator["singlet_triplet"]
    doublet = operator["doublet"]
    kramers = operator["kramers"]
    zero_maximum = max(
        value
        for family in operator["zero_SOC_equivalence"].values()
        for value in family.values()
    )
    st_states = SingletTripletSOCConfigV220().model_space().states
    doublet_states = DoubletSOCConfigV220().model_space().states
    st_projectors = AnalyticSingletTripletSOCProviderV220().projectors
    doublet_projectors = AnalyticDoubletSOCProviderV220().projectors
    st_exact = exact["singlet_triplet"]
    doublet_exact = exact["doublet"]
    maximum_gaussian_grid = max(
        row["absolute_error"]
        for row in exact["gaussian_exact_grid_population_errors"].values()
    )
    maximum_dense_restart = max(restart["dense_errors"].values())
    maximum_sparse_restart = max(restart["sparse_errors"].values())
    maximum_gauge_dynamics = max(gauge["dynamics_errors"].values())

    new_checks = {
        # Shared physical-SOC contract: 8
        "shared::H_decomposition_both_models": max(
            st["H_composition_error"], doublet["H_composition_error"]
        ) <= thresholds.maximum_composition_error,
        "shared::K_decomposition_both_models": max(
            st["K_composition_error"], doublet["K_composition_error"]
        ) <= thresholds.maximum_composition_error,
        "shared::structure_and_time_reversal": max(
            st["maximum_time_reversal_residual"],
            doublet["maximum_time_reversal_residual"],
        ) <= thresholds.maximum_time_reversal_residual,
        "shared::cross_geometry_differentials": (
            st["differential_report"]["passed"]
            and doublet["differential_report"]["passed"]
        ),
        "shared::zero_SOC_operator_equivalence": zero_maximum == 0.0,
        "shared::moving_complex_gauge_differentials": gauge["differential"]["passed"],
        "shared::SOC_force_derivatives": max(
            st["soc_force_error"], doublet["soc_force_error"]
        ) <= thresholds.maximum_force_error,
        "shared::wrong_SOC_derivative_detected": (
            not operator["wrong_SOC_derivative"]["checks"]["K_decomposition"]
            and not operator["wrong_SOC_derivative"]["checks"]["SOC_force_derivative"]
        ),
        # Singlet-triplet model: 6
        "singlet_triplet::complete_multiplets": (
            len(st_states) == 4
            and [state.multiplicity for state in st_states] == [1, 3, 3, 3]
        ),
        "singlet_triplet::even_time_reversal_square": st[
            "time_reversal_square_residual"
        ] <= thresholds.maximum_time_reversal_residual,
        "singlet_triplet::constant_SOC_zero_derivative": (
            operator["constant_SOC"]["H_soc_norm"] > 0.0
            and operator["constant_SOC"]["K_soc_norm"] == 0.0
        ),
        "singlet_triplet::coordinate_SOC_H_and_K_nonzero": (
            np.linalg.norm(AnalyticSingletTripletSOCProviderV220().components([0.17]).H_soc) > 0.0
            and np.linalg.norm(AnalyticSingletTripletSOCProviderV220().components([0.17]).K_soc) > 0.0
        ),
        "singlet_triplet::projector_completeness": np.array_equal(
            sum(st_projectors.values()), np.eye(4)
        ),
        "singlet_triplet::physical_population_transfer": st_exact[
            "final_populations"
        ]["triplet"] > 5.0e-5,
        # Doublet/Kramers model: 8
        "doublet::two_complete_doublets": (
            len(doublet_states) == 4
            and [state.source_root for state in doublet_states].count("D1") == 2
            and [state.source_root for state in doublet_states].count("D2") == 2
        ),
        "doublet::fermionic_time_reversal_square": kramers[
            "time_reversal_square_residual"
        ] <= thresholds.maximum_time_reversal_residual,
        "doublet::H_time_reversal": kramers[
            "maximum_H_time_reversal_residual"
        ] <= thresholds.maximum_time_reversal_residual,
        "doublet::K_time_reversal": kramers[
            "maximum_K_time_reversal_residual"
        ] <= thresholds.maximum_time_reversal_residual,
        "doublet::kramers_pair_degeneracy": kramers["maximum_pair_splitting"]
        <= thresholds.maximum_kramers_splitting,
        "doublet::transformed_time_reversal_operator": gauge[
            "transformed_time_reversal_residual"
        ] <= thresholds.maximum_time_reversal_residual,
        "doublet::projector_gauge_invariance": gauge[
            "maximum_projector_population_error"
        ] <= thresholds.maximum_gauge_error,
        "doublet::broken_kramers_fixture_detected": (
            operator["broken_kramers"]["time_reversal_residual"] > 1.0e-4
            and operator["broken_kramers"]["pair_splitting"] > 1.0e-5
            and np.array_equal(sum(doublet_projectors.values()), np.eye(4))
        ),
        # Independent exact-grid and convergence: 5
        "exact_grid::norm_conservation_both_models": max(
            st_exact["maximum_norm_drift"], doublet_exact["maximum_norm_drift"]
        ) <= thresholds.maximum_grid_norm_drift,
        "exact_grid::energy_conservation_both_models": max(
            st_exact["maximum_energy_drift"], doublet_exact["maximum_energy_drift"]
        ) <= thresholds.maximum_grid_energy_drift,
        "exact_grid::second_order_timestep_convergence": exact[
            "timestep_convergence"
        ]["observed_order"] >= thresholds.minimum_timestep_order,
        "exact_grid::grid_spacing_and_box_convergence": exact[
            "grid_resolution_and_box_population_error"
        ] <= thresholds.maximum_grid_population_resolution_error,
        "exact_grid::short_time_gaussian_population_agreement": maximum_gaussian_grid
        <= thresholds.maximum_gaussian_grid_population_error,
        # SOC restart and failure controls: 5
        "restart::SOC_dense_segment_equivalence": maximum_dense_restart
        <= thresholds.maximum_restart_error,
        "restart::SOC_sparse_segment_equivalence": (
            maximum_sparse_restart <= thresholds.maximum_restart_error
            and restart["sparse_edges"]["checkpoint"] == [[3, 8]]
            and restart["sparse_edges"]["uninterrupted"]
            == restart["sparse_edges"]["resumed"]
        ),
        "restart::moving_complex_gauge_dynamics": maximum_gauge_dynamics
        <= thresholds.maximum_gauge_error,
        "restart::changed_SOC_provenance_rejected": restart[
            "changed_soc_provenance_rejected"
        ],
        "restart::SOC_checkpoint_corruption_rejected": restart[
            "corruption_rejected"
        ],
    }
    if len(new_checks) != 32:
        raise AssertionError("v0.22.0 campaign must define exactly 32 new gates.")
    new_checks = {name: bool(passed) for name, passed in new_checks.items()}
    inherited_checks = {
        f"inherited_v0214::{name}": bool(passed)
        for name, passed in inherited["acceptance"]["checks"].items()
    }
    checks = {**inherited_checks, **new_checks}
    if len(checks) != 53:
        raise AssertionError("v0.22.0 campaign must define exactly 53 total gates.")
    return {
        "release": "v0.22.0",
        "theme": "first physical analytic SOC with singlet-triplet and Kramers-doublet references",
        "operator_contract": operator,
        "gauge_covariance": gauge,
        "exact_grid": exact,
        "checkpoint_restart": restart,
        "inherited_v0214": inherited,
        "soc": {
            "physical_hamiltonian_introduced": True,
            "physical_derivative_introduced": True,
            "analytic_models_only": True,
            "singlet_triplet_model": True,
            "kramers_doublet_model": True,
            "spin_free_mode_permanent": True,
            "ab_initio_SOC_validated": False,
            "external_magnetic_field": False,
        },
        "pyscf": {
            "installed_in_build_environment": bool(importlib.util.find_spec("pyscf")),
            "runtime_validated": False,
            "note": "v0.22.0 validates physical analytic SOC models only.",
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
