"""Canonical v0.21.4 pre-SOC differential and restart certification campaign."""

from dataclasses import asdict, dataclass
import importlib.util
from pathlib import Path
import tempfile
import numpy as np

from .block_sparse_molecular_v21 import BlockMolecularTBFV21, BlockSparseSettingsV21
from .checkpoint_restart_v214 import (
    SelfConsistentBlockCheckpointV214,
    SelfConsistentBlockSettingsV214,
    load_self_consistent_checkpoint_v214,
    run_self_consistent_block_dynamics_v214,
    save_self_consistent_checkpoint_v214,
    settings_fingerprint_v214,
)
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from .density_guidance_v213 import BlockDensityMatrixGuidanceV213
from .electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from .electronic_operator_v21 import ElectronicOperatorPointV21, ElectronicOperatorSnapshotV21
from .provider_differential_audit_v214 import audit_provider_differentials_v214
from .synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)
from .v213_benchmark import run_v0213_release_benchmark
from .zero_soc_rehearsal_v214 import (
    ZeroSOCRehearsalProviderV214,
    audit_zero_soc_equivalence_v214,
)


@dataclass(frozen=True)
class V214AcceptanceThresholds:
    max_fixed_K_differential_error: float = 1.0e-10
    max_fixed_D_differential_error: float = 1.0e-14
    max_gauge_K_differential_error: float = 2.0e-10
    max_gauge_D_differential_error: float = 2.0e-9
    max_overlap_isometry_residual: float = 1.0e-12
    max_zero_soc_operator_error: float = 0.0
    max_zero_soc_dynamics_error: float = 1.0e-13
    max_restart_position_error: float = 1.0e-13
    max_restart_momentum_error: float = 1.0e-13
    max_restart_coefficient_error: float = 2.0e-12
    max_moving_frame_restart_error: float = 2.0e-12
    require_inherited_v0213: bool = True


def _model_space_v214(representation):
    return ElectronicModelSpaceV213(
        name=f"v0.21.4 two-state {representation} certification space",
        representation=representation,
        states=(
            ElectronicStateDescriptorV213("state-0"),
            ElectronicStateDescriptorV213("state-1"),
        ),
    ).validate()


def _provenance_v214(representation="fixed_general"):
    return ElectronicOperatorProvenanceV213(
        model_name="v0.21.4 zero-SOC integration rehearsal",
        model_version="1",
        model_space=_model_space_v214(representation),
        spin_free_method="analytic complex linear fixture",
        soc_enabled=False,
        soc_method="none",
        derivative_method="analytic physical operator derivative",
        parameters={
            "base_seed": 21450,
            "gauge_seed": 21451,
            "phase_gradient": [[0.21], [-0.16]],
        },
    ).validate()


def _base_provider_v214():
    return SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2,
            nq=1,
            mass=28.0,
            seed=21450,
            base_scale=0.025,
            derivative_scale=0.008,
        )
    )


def _contracted_provider_v214(provenance):
    return ContractedElectronicOperatorProviderV213(
        _base_provider_v214(), provenance
    )


def _gauge_v214():
    return PhaseMixingGaugeV21(
        random_unitary_v21(2, 21451),
        np.asarray([[0.21], [-0.16]]),
        np.asarray([0.18, -0.27]),
    )


def _basis_v214():
    return [
        BlockMolecularTBFV21(
            3, np.asarray([-0.55]), np.asarray([0.22]), np.asarray([[1.15]])
        ),
        BlockMolecularTBFV21(
            8, np.asarray([0.45]), np.asarray([-0.12]), np.asarray([[1.45]])
        ),
    ]


def _coefficients_v214():
    return np.asarray([0.72 + 0.10j, 0.18 - 0.24j, -0.11 + 0.28j, 0.31 + 0.05j])


def _phase_aligned_metric_error(reference, candidate, metric):
    overlap = np.vdot(reference, metric @ candidate)
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else np.exp(-1.0j * np.angle(overlap))
    difference = phase * candidate - reference
    return float(np.sqrt(max(np.real(np.vdot(difference, metric @ difference)), 0.0)))


def _trajectory_errors(reference, candidate):
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
        "coefficient": _phase_aligned_metric_error(
            reference["final_coefficients"],
            candidate["final_coefficients"],
            reference["final_S"],
        ),
    }


class _PerturbedProviderV214:
    def __init__(self, base, *, perturb_K=False, erase_D=False):
        self.base = base
        self.perturb_K = bool(perturb_K)
        self.erase_D = bool(erase_D)

    def evaluate_snapshot(self, q):
        snapshot = self.base.evaluate_snapshot(q)
        K = snapshot.point.dH_dq.copy()
        D = snapshot.point.connection_q.copy()
        if self.perturb_K:
            K[0] += 1.0e-3 * np.eye(snapshot.point.nstate)
        if self.erase_D:
            D[:] = 0.0
        point = ElectronicOperatorPointV21(
            q=snapshot.point.q.copy(),
            H=snapshot.point.H.copy(),
            dH_dq=K,
            connection_q=D,
            mass_matrix_q_au=snapshot.point.mass_matrix_q_au.copy(),
            metadata=dict(snapshot.point.metadata),
        ).validate()
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=snapshot.state_vectors.copy(),
            parent_snapshot=snapshot,
        ).validate()

    def snapshot_overlap(self, left, right):
        return self.base.snapshot_overlap(left.parent_snapshot, right.parent_snapshot)


def _differential_campaign_v214():
    q = np.asarray([0.17])
    fixed_provenance = _provenance_v214("fixed_general")
    fixed = ContractedElectronicOperatorProviderV213(
        _base_provider_v214(), fixed_provenance
    )
    fixed_report = audit_provider_differentials_v214(fixed, q, fixed_provenance)

    local_provenance = _provenance_v214("local_general")
    local = ContractedElectronicOperatorProviderV213(
        GaugeTransformedOperatorProviderV21(_base_provider_v214(), _gauge_v214()),
        local_provenance,
    )
    gauge_report = audit_provider_differentials_v214(local, q, local_provenance)
    wrong_K = audit_provider_differentials_v214(
        _PerturbedProviderV214(fixed, perturb_K=True), q, fixed_provenance
    )
    wrong_D = audit_provider_differentials_v214(
        _PerturbedProviderV214(local, erase_D=True), q, local_provenance
    )
    return {
        "fixed_frame": fixed_report.as_dict(),
        "coordinate_dependent_complex_frame": gauge_report.as_dict(),
        "wrong_K_fixture": wrong_K.as_dict(),
        "wrong_D_fixture": wrong_D.as_dict(),
    }


def _zero_soc_campaign_v214():
    provenance = _provenance_v214("fixed_general")
    rehearsal = ZeroSOCRehearsalProviderV214(_base_provider_v214(), provenance)
    equivalence = audit_zero_soc_equivalence_v214(
        _base_provider_v214(),
        rehearsal,
        [np.asarray([-0.3]), np.asarray([0.0]), np.asarray([0.4])],
        tolerance=0.0,
    )
    differential = audit_provider_differentials_v214(
        rehearsal, np.asarray([0.17]), provenance
    )
    settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    common = dict(
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        dt=0.002,
        steps=6,
        store_every=2,
        settings=settings,
    )
    spin_free = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance), provenance, **common
    )
    zero_soc = run_self_consistent_block_dynamics_v214(
        ZeroSOCRehearsalProviderV214(_base_provider_v214(), provenance),
        provenance,
        **common,
    )
    return {
        "operator_equivalence": equivalence.as_dict(),
        "differential_contract": differential.as_dict(),
        "dynamics_errors": _trajectory_errors(spin_free, zero_soc),
    }


def _restart_campaign_v214():
    provenance = _provenance_v214("fixed_general")
    dense_settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    common = dict(dt=0.002, store_every=2, settings=dense_settings)
    full = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        steps=10,
        **common,
    )
    first = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        steps=4,
        **common,
    )
    with tempfile.TemporaryDirectory(prefix="v214-checkpoint-") as directory:
        path = save_self_consistent_checkpoint_v214(
            Path(directory) / "trajectory.npz", first["checkpoint"]
        )
        loaded = load_self_consistent_checkpoint_v214(
            path,
            expected_provider_fingerprint=provenance.fingerprint(),
            expected_settings_fingerprint=settings_fingerprint_v214(dense_settings),
        )
        with np.load(path, allow_pickle=False) as archive:
            corrupted_arrays = {name: archive[name].copy() for name in archive.files}
        corrupted_arrays["q"][0, 0] += 1.0e-3
        corrupted_path = Path(directory) / "corrupted.npz"
        with corrupted_path.open("wb") as handle:
            np.savez_compressed(handle, **corrupted_arrays)
        corruption_rejected = False
        try:
            load_self_consistent_checkpoint_v214(corrupted_path)
        except ValueError:
            corruption_rejected = True
    resumed = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        checkpoint=loaded,
        steps=6,
        store_every=2,
        settings=dense_settings,
    )

    sparse_graph = BlockSparseSettingsV21(
        enter_score=0.30,
        exit_score=0.10,
        search_overlap_floor=1.0e-8,
        local_omitted_score_l2_budget=1.0,
        use_kdtree=False,
    )
    sparse_settings = SelfConsistentBlockSettingsV214(
        graph=sparse_graph,
        use_dense_reference=False,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    sparse_full = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        dt=0.002,
        steps=10,
        store_every=2,
        settings=sparse_settings,
    )
    sparse_first = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        dt=0.002,
        steps=4,
        store_every=2,
        settings=sparse_settings,
    )
    sparse_resumed = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        checkpoint=sparse_first["checkpoint"],
        steps=6,
        store_every=2,
        settings=sparse_settings,
    )

    local_provenance = _provenance_v214("local_general")
    gauge = _gauge_v214()

    def gauge_provider():
        return ContractedElectronicOperatorProviderV213(
            GaugeTransformedOperatorProviderV21(_base_provider_v214(), gauge),
            local_provenance,
        )

    transformed_C = np.concatenate(
        [
            gauge.matrix(item.q).conj().T @ _coefficients_v214()[2 * i : 2 * i + 2]
            for i, item in enumerate(_basis_v214())
        ]
    )
    moving_full = run_self_consistent_block_dynamics_v214(
        gauge_provider(),
        local_provenance,
        initial_basis=_basis_v214(),
        C0=transformed_C,
        dt=0.002,
        steps=8,
        store_every=2,
        settings=dense_settings,
    )
    moving_first = run_self_consistent_block_dynamics_v214(
        gauge_provider(),
        local_provenance,
        initial_basis=_basis_v214(),
        C0=transformed_C,
        dt=0.002,
        steps=3,
        store_every=1,
        settings=dense_settings,
    )
    moving_resumed = run_self_consistent_block_dynamics_v214(
        gauge_provider(),
        local_provenance,
        checkpoint=moving_first["checkpoint"],
        steps=5,
        store_every=1,
        settings=dense_settings,
    )

    guidance = BlockDensityMatrixGuidanceV213(dense_settings.guidance)
    guide_basis = _basis_v214()
    guide_provider = _contracted_provider_v214(provenance)
    guidance.on_insert(
        guide_basis[0], guide_provider, guide_density=np.diag([1.0, 0.0])
    )
    guidance.on_insert(
        guide_basis[1], guide_provider, guide_density=np.diag([0.0, 1.0])
    )
    guide_checkpoint = SelfConsistentBlockCheckpointV214.create(
        step=3,
        dt=0.002,
        provider_fingerprint=provenance.fingerprint(),
        settings_fingerprint=settings_fingerprint_v214(dense_settings),
        basis=guide_basis,
        coefficients=np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex),
        nstate=2,
        guidance=guidance,
    )
    guide_resumed = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        checkpoint=guide_checkpoint,
        steps=0,
        store_every=1,
        settings=dense_settings,
    )
    return {
        "dense_errors": _trajectory_errors(full, resumed),
        "checkpoint_roundtrip": {
            "step": loaded.step,
            "time": loaded.time,
            "integrity_digest": loaded.integrity_digest,
            "digest_recomputed_exactly": (
                loaded.integrity_digest == loaded.computed_integrity_digest()
            ),
            "corruption_rejected": corruption_rejected,
        },
        "sparse_errors": _trajectory_errors(sparse_full, sparse_resumed),
        "sparse_edges": {
            "checkpoint": sparse_first["checkpoint"].active_uid_edges.tolist(),
            "uninterrupted_final": [list(edge) for edge in sparse_full["final_active_uid_edges"]],
            "resumed_final": [list(edge) for edge in sparse_resumed["final_active_uid_edges"]],
        },
        "moving_frame_errors": _trajectory_errors(moving_full, moving_resumed),
        "retained_zero_block_density": {
            "guide_mask": guide_resumed["checkpoint"].guide_mask.tolist(),
            "second_density": guide_resumed["checkpoint"].guide_densities[1].tolist(),
            "retained_density_uses": guide_resumed["guidance_diagnostics"][
                "retained_density_uses"
            ],
        },
    }


def _lifecycle_campaign_v214():
    provenance = _provenance_v214("fixed_general")
    settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    child = BlockMolecularTBFV21(
        13, np.asarray([1.15]), np.asarray([0.04]), np.asarray([[1.25]])
    )

    def adaptation(global_step, basis, coefficients, metric):
        if global_step == 2:
            return {"insert": child, "guide_parent_uid": 8}
        if global_step == 3:
            return {"prune_index": 2}
        return None

    first = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        initial_basis=_basis_v214(),
        C0=_coefficients_v214(),
        dt=0.002,
        steps=2,
        store_every=1,
        settings=settings,
        adaptation_policy=adaptation,
    )
    resumed = run_self_consistent_block_dynamics_v214(
        _contracted_provider_v214(provenance),
        provenance,
        checkpoint=first["checkpoint"],
        steps=2,
        store_every=1,
        settings=settings,
        adaptation_policy=adaptation,
    )
    return {
        "checkpoint_uids": first["checkpoint"].uids.tolist(),
        "checkpoint_guide_mask": first["checkpoint"].guide_mask.tolist(),
        "resumed_events": [
            {"step": event["step"], "kind": event["kind"]}
            for event in resumed["adaptation_events"]
        ],
        "final_uids": [item.uid for item in resumed["final_basis"]],
    }


def run_v0214_release_benchmark():
    differential = _differential_campaign_v214()
    zero_soc = _zero_soc_campaign_v214()
    restart = _restart_campaign_v214()
    lifecycle = _lifecycle_campaign_v214()
    inherited = bool(run_v0213_release_benchmark()["acceptance"]["passed"])
    thresholds = V214AcceptanceThresholds()
    fixed = differential["fixed_frame"]
    gauge = differential["coordinate_dependent_complex_frame"]
    zero_operator = zero_soc["operator_equivalence"]
    zero_dynamics = zero_soc["dynamics_errors"]
    dense = restart["dense_errors"]
    sparse = restart["sparse_errors"]
    moving = restart["moving_frame_errors"]
    checks = {
        "fixed_frame_differential_contract": (
            fixed["passed"]
            and fixed["maximum_hamiltonian_derivative_scaled_error"]
            <= thresholds.max_fixed_K_differential_error
            and fixed["maximum_connection_scaled_error"]
            <= thresholds.max_fixed_D_differential_error
        ),
        "complex_gauge_differential_contract": (
            gauge["passed"]
            and gauge["maximum_hamiltonian_derivative_scaled_error"]
            <= thresholds.max_gauge_K_differential_error
            and gauge["maximum_connection_scaled_error"]
            <= thresholds.max_gauge_D_differential_error
            and gauge["maximum_overlap_isometry_residual"]
            <= thresholds.max_overlap_isometry_residual
        ),
        "wrong_K_detected": (
            differential["wrong_K_fixture"]["checks"]["structural_invariants"]
            and not differential["wrong_K_fixture"]["checks"][
                "physical_H_derivatives"
            ]
        ),
        "wrong_D_detected": (
            differential["wrong_D_fixture"]["checks"]["structural_invariants"]
            and not differential["wrong_D_fixture"]["checks"][
                "derivative_connections"
            ]
        ),
        "zero_soc_H_equivalence": zero_operator["maximum_H_error"]
        <= thresholds.max_zero_soc_operator_error,
        "zero_soc_K_equivalence": zero_operator["maximum_K_error"]
        <= thresholds.max_zero_soc_operator_error,
        "zero_soc_D_mass_overlap_equivalence": max(
            zero_operator["maximum_D_error"],
            zero_operator["maximum_mass_error"],
            zero_operator["maximum_overlap_error"],
        )
        <= thresholds.max_zero_soc_operator_error,
        "zero_soc_differential_contract": zero_soc["differential_contract"][
            "passed"
        ],
        "zero_soc_dynamics_positions": zero_dynamics["position"]
        <= thresholds.max_zero_soc_dynamics_error,
        "zero_soc_dynamics_momenta": zero_dynamics["momentum"]
        <= thresholds.max_zero_soc_dynamics_error,
        "zero_soc_dynamics_coefficients": zero_dynamics["coefficient"]
        <= thresholds.max_zero_soc_dynamics_error,
        "dense_restart_positions": dense["position"]
        <= thresholds.max_restart_position_error,
        "dense_restart_momenta": dense["momentum"]
        <= thresholds.max_restart_momentum_error,
        "dense_restart_coefficients": dense["coefficient"]
        <= thresholds.max_restart_coefficient_error,
        "checkpoint_integrity_roundtrip": restart["checkpoint_roundtrip"][
            "digest_recomputed_exactly"
        ],
        "checkpoint_corruption_rejected": restart["checkpoint_roundtrip"][
            "corruption_rejected"
        ],
        "sparse_graph_restart": (
            sparse["position"] <= thresholds.max_restart_position_error
            and sparse["momentum"] <= thresholds.max_restart_momentum_error
            and sparse["coefficient"] <= thresholds.max_restart_coefficient_error
            and restart["sparse_edges"]["checkpoint"] == [[3, 8]]
            and restart["sparse_edges"]["uninterrupted_final"]
            == restart["sparse_edges"]["resumed_final"]
        ),
        "moving_complex_frame_restart": max(moving.values())
        <= thresholds.max_moving_frame_restart_error,
        "zero_block_guide_density_restart": (
            restart["retained_zero_block_density"]["guide_mask"] == [True, True]
            and restart["retained_zero_block_density"]["retained_density_uses"] >= 1
        ),
        "global_adaptive_lifecycle_restart": (
            lifecycle["checkpoint_uids"] == [3, 8, 13]
            and lifecycle["checkpoint_guide_mask"] == [True, True, True]
            and lifecycle["resumed_events"] == [{"step": 3, "kind": "prune"}]
            and lifecycle["final_uids"] == [3, 8]
        ),
        "inherited_v0213": inherited if thresholds.require_inherited_v0213 else True,
    }
    return {
        "release": "v0.21.4",
        "theme": "pre-SOC differential provider and deterministic restart certification",
        "provider_differential_certification": differential,
        "explicit_zero_soc_rehearsal": zero_soc,
        "checkpoint_restart": restart,
        "adaptive_lifecycle_restart": lifecycle,
        "inherited_v0213_acceptance": inherited,
        "soc": {
            "physical_hamiltonian_introduced": False,
            "physical_derivative_introduced": False,
            "spin_free_mode_permanent": True,
            "first_physical_soc_target": "v0.22",
        },
        "pyscf": {
            "installed_in_build_environment": bool(
                importlib.util.find_spec("pyscf") is not None
            ),
            "runtime_validated": False,
            "note": "v0.21.4 certifies analytic interfaces and restart behavior only.",
        },
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "thresholds": asdict(thresholds),
        },
    }
