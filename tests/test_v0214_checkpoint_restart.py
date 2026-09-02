from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.block_sparse_molecular_v21 import (
    BlockMolecularTBFV21,
    BlockSparseMolecularGraphV21,
    BlockSparseSettingsV21,
)
from gaussian_dynamics.checkpoint_restart_v214 import (
    SelfConsistentBlockCheckpointV214,
    SelfConsistentBlockSettingsV214,
    load_self_consistent_checkpoint_v214,
    run_self_consistent_block_dynamics_v214,
    save_self_consistent_checkpoint_v214,
    settings_fingerprint_v214,
)
from gaussian_dynamics.density_guidance_v213 import BlockDensityMatrixGuidanceV213
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ContractedElectronicOperatorProviderV213,
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)


def _provenance(parameter=1.0):
    return ElectronicOperatorProvenanceV213(
        model_name="v0.21.4 restart fixture",
        model_version="1",
        model_space=ElectronicModelSpaceV213(
            name="two-state fixed restart space",
            representation="fixed_general",
            states=(
                ElectronicStateDescriptorV213("state-0"),
                ElectronicStateDescriptorV213("state-1"),
            ),
        ),
        spin_free_method="analytic linear fixture",
        parameters={"parameter": parameter},
    )


def _provider(provenance=None):
    provenance = _provenance() if provenance is None else provenance
    base = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2,
            nq=1,
            mass=28.0,
            seed=21420,
            base_scale=0.025,
            derivative_scale=0.008,
        )
    )
    return ContractedElectronicOperatorProviderV213(base, provenance)


def _basis():
    return [
        BlockMolecularTBFV21(
            3, np.asarray([-0.55]), np.asarray([0.22]), np.asarray([[1.15]])
        ),
        BlockMolecularTBFV21(
            8, np.asarray([0.45]), np.asarray([-0.12]), np.asarray([[1.45]])
        ),
    ]


def _coefficients():
    return np.asarray([0.72 + 0.10j, 0.18 - 0.24j, -0.11 + 0.28j, 0.31 + 0.05j])


def _phase_aligned_metric_error(reference, candidate, metric):
    overlap = np.vdot(reference, metric @ candidate)
    phase = 1.0 + 0.0j if abs(overlap) < 1.0e-30 else np.exp(-1.0j * np.angle(overlap))
    difference = phase * candidate - reference
    return float(
        np.sqrt(max(np.real(np.vdot(difference, metric @ difference)), 0.0))
    )


def test_dense_checkpoint_restart_matches_uninterrupted_trajectory(tmp_path):
    provenance = _provenance()
    settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    dt = 0.002
    full = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=dt,
        steps=10,
        store_every=2,
        settings=settings,
    )
    first = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=dt,
        steps=4,
        store_every=2,
        settings=settings,
    )
    path = save_self_consistent_checkpoint_v214(
        tmp_path / "restart-v214.npz", first["checkpoint"]
    )
    loaded = load_self_consistent_checkpoint_v214(
        path,
        expected_provider_fingerprint=provenance.fingerprint(),
        expected_settings_fingerprint=settings_fingerprint_v214(settings),
    )
    resumed = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        checkpoint=loaded,
        steps=6,
        store_every=2,
        settings=settings,
    )

    q_error = max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    )
    p_error = max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    )
    coefficient_error = _phase_aligned_metric_error(
        full["final_coefficients"], resumed["final_coefficients"], full["final_S"]
    )
    assert q_error < 1.0e-13
    assert p_error < 1.0e-13
    assert coefficient_error < 2.0e-12
    assert loaded.integrity_digest == loaded.computed_integrity_digest()
    assert resumed["checkpoint"].step == 10
    assert resumed["checkpoint"].time == 10 * dt
    assert resumed["records"][0]["step"] == 4
    assert resumed["records"][-1]["step"] == 10
    assert resumed["restart_source"] is True


def test_checkpoint_rejects_provenance_settings_and_integrity_mismatches(tmp_path):
    provenance = _provenance()
    settings = SelfConsistentBlockSettingsV214()
    run = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=1,
        store_every=1,
        settings=settings,
    )
    path = save_self_consistent_checkpoint_v214(
        tmp_path / "identity-v214.npz", run["checkpoint"]
    )
    with pytest.raises(ValueError, match="provider provenance"):
        load_self_consistent_checkpoint_v214(
            path, expected_provider_fingerprint=_provenance(2.0).fingerprint()
        )
    changed_settings = SelfConsistentBlockSettingsV214(momentum_tolerance=2.0e-10)
    with pytest.raises(ValueError, match="propagation-settings"):
        run_self_consistent_block_dynamics_v214(
            _provider(provenance),
            provenance,
            checkpoint=run["checkpoint"],
            steps=1,
            settings=changed_settings,
        )

    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name].copy() for name in archive.files}
    arrays["q"][0, 0] += 1.0e-3
    corrupted = tmp_path / "corrupted-v214.npz"
    with corrupted.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    with pytest.raises(ValueError, match="integrity digest"):
        load_self_consistent_checkpoint_v214(corrupted)

    reversed_edge = replace(
        run["checkpoint"],
        active_uid_edges=np.asarray([[8, 3]], dtype=np.int64),
        integrity_digest="0" * 64,
    )
    reversed_edge = replace(
        reversed_edge,
        integrity_digest=reversed_edge.computed_integrity_digest(),
    )
    with pytest.raises(ValueError, match="ascending uid order"):
        reversed_edge.validate()

    with np.load(path, allow_pickle=False) as archive:
        byte_manifest_arrays = {
            name: archive[name].copy() for name in archive.files
        }
    byte_manifest_arrays["manifest_json"] = np.asarray(
        byte_manifest_arrays["manifest_json"].item().encode("utf-8")
    )
    byte_manifest = tmp_path / "byte-manifest-v214.npz"
    with byte_manifest.open("wb") as handle:
        np.savez_compressed(handle, **byte_manifest_arrays)
    assert load_self_consistent_checkpoint_v214(byte_manifest).step == 1


def test_v0214_runner_requires_provider_emitted_provenance_identity():
    provenance = _provenance()
    bare_provider = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(
            nstate=2,
            nq=1,
            mass=28.0,
            seed=21420,
            base_scale=0.025,
            derivative_scale=0.008,
        )
    )

    with pytest.raises(ValueError, match="emit the declared provenance fingerprint"):
        run_self_consistent_block_dynamics_v214(
            bare_provider,
            provenance,
            initial_basis=_basis(),
            C0=_coefficients(),
            dt=0.002,
            steps=0,
        )


def test_checkpoint_preserves_guide_density_for_an_exact_zero_local_block():
    provenance = _provenance()
    provider = _provider(provenance)
    settings = SelfConsistentBlockSettingsV214()
    basis = _basis()
    guidance = BlockDensityMatrixGuidanceV213(settings.guidance)
    guidance.on_insert(
        basis[0], provider, guide_density=np.diag([1.0, 0.0]).astype(complex)
    )
    guidance.on_insert(
        basis[1], provider, guide_density=np.diag([0.0, 1.0]).astype(complex)
    )
    checkpoint = SelfConsistentBlockCheckpointV214.create(
        step=3,
        dt=0.002,
        provider_fingerprint=provenance.fingerprint(),
        settings_fingerprint=settings_fingerprint_v214(settings),
        basis=basis,
        coefficients=np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex),
        nstate=2,
        guidance=guidance,
    )
    resumed = run_self_consistent_block_dynamics_v214(
        provider,
        provenance,
        checkpoint=checkpoint,
        steps=0,
        store_every=1,
        settings=settings,
    )

    assert np.array_equal(resumed["checkpoint"].guide_mask, [True, True])
    assert np.allclose(
        resumed["checkpoint"].guide_densities[1], np.diag([0.0, 1.0])
    )
    assert resumed["guidance_diagnostics"]["retained_density_uses"] >= 1


def test_sparse_graph_hysteresis_state_has_a_validated_public_restart_path():
    provider = _provider()
    settings = BlockSparseSettingsV21(
        enter_score=0.05,
        exit_score=0.02,
        search_overlap_floor=1.0e-8,
        use_kdtree=False,
    )
    graph = BlockSparseMolecularGraphV21(provider, 0.002, settings)
    assert graph.restore_active_uid_edges_v214([(8, 3)], [3, 8]) == ((3, 8),)
    assert graph.active_uid_edges_v214() == ((3, 8),)
    with pytest.raises(ValueError, match="non-live"):
        graph.restore_active_uid_edges_v214([(3, 99)], [3, 8])
    with pytest.raises(ValueError, match="unique"):
        graph.restore_active_uid_edges_v214([(3, 8), (8, 3)], [3, 8])


def test_sparse_checkpoint_restart_restores_hysteretic_uid_edges_and_trajectory():
    provenance = _provenance()
    graph = BlockSparseSettingsV21(
        enter_score=0.30,
        exit_score=0.10,
        search_overlap_floor=1.0e-8,
        local_omitted_score_l2_budget=1.0,
        use_kdtree=False,
    )
    settings = SelfConsistentBlockSettingsV214(
        graph=graph,
        use_dense_reference=False,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    full = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=10,
        store_every=2,
        settings=settings,
    )
    first = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=4,
        store_every=2,
        settings=settings,
    )
    resumed = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        checkpoint=first["checkpoint"],
        steps=6,
        store_every=2,
        settings=settings,
    )

    assert first["checkpoint"].active_uid_edges.tolist() == [[3, 8]]
    assert full["final_active_uid_edges"] == resumed["final_active_uid_edges"] == (
        (3, 8),
    )
    assert max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) < 1.0e-13
    assert max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) < 1.0e-13
    assert _phase_aligned_metric_error(
        full["final_coefficients"], resumed["final_coefficients"], full["final_S"]
    ) < 2.0e-12


def test_restart_uses_global_adaptation_steps_and_preserves_block_lifecycle():
    provenance = _provenance()
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

    full = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=4,
        store_every=1,
        settings=settings,
        adaptation_policy=adaptation,
    )
    first = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        initial_basis=_basis(),
        C0=_coefficients(),
        dt=0.002,
        steps=2,
        store_every=1,
        settings=settings,
        adaptation_policy=adaptation,
    )
    resumed = run_self_consistent_block_dynamics_v214(
        _provider(provenance),
        provenance,
        checkpoint=first["checkpoint"],
        steps=2,
        store_every=1,
        settings=settings,
        adaptation_policy=adaptation,
    )

    assert first["checkpoint"].uids.tolist() == [3, 8, 13]
    assert first["checkpoint"].guide_mask.tolist() == [True, True, True]
    assert [(event["step"], event["kind"]) for event in full["adaptation_events"]] == [
        (2, "insert"),
        (3, "prune"),
    ]
    assert [(event["step"], event["kind"]) for event in resumed["adaptation_events"]] == [
        (3, "prune")
    ]
    assert [item.uid for item in full["final_basis"]] == [3, 8]
    assert [item.uid for item in resumed["final_basis"]] == [3, 8]
    assert _phase_aligned_metric_error(
        full["final_coefficients"], resumed["final_coefficients"], full["final_S"]
    ) < 2.0e-12


def test_checkpoint_restart_is_stable_in_a_coordinate_dependent_complex_frame():
    provenance = ElectronicOperatorProvenanceV213(
        model_name="v0.21.4 moving-frame restart fixture",
        model_version="1",
        model_space=ElectronicModelSpaceV213(
            name="two-state local restart space",
            representation="local_general",
            states=(
                ElectronicStateDescriptorV213("state-0"),
                ElectronicStateDescriptorV213("state-1"),
            ),
        ),
        spin_free_method="analytic linear fixture",
        parameters={
            "base_seed": 21420,
            "gauge_seed": 21441,
            "phase_gradient": [[0.21], [-0.16]],
        },
    )
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(2, 21441),
        np.asarray([[0.21], [-0.16]]),
        np.asarray([0.18, -0.27]),
    )

    def provider():
        return ContractedElectronicOperatorProviderV213(
            GaugeTransformedOperatorProviderV21(
                SyntheticLinearOperatorProviderV21(
                    SyntheticLinearOperatorConfigV21(
                        nstate=2,
                        nq=1,
                        mass=28.0,
                        seed=21420,
                        base_scale=0.025,
                        derivative_scale=0.008,
                    )
                ),
                gauge,
            ),
            provenance,
        )

    initial_basis = _basis()
    base_coefficients = _coefficients()
    transformed_coefficients = np.concatenate(
        [
            gauge.matrix(item.q).conj().T @ base_coefficients[2 * index : 2 * index + 2]
            for index, item in enumerate(initial_basis)
        ]
    )
    settings = SelfConsistentBlockSettingsV214(
        use_dense_reference=True,
        corrector_iterations=3,
        momentum_tolerance=1.0e-12,
    )
    full = run_self_consistent_block_dynamics_v214(
        provider(),
        provenance,
        initial_basis=initial_basis,
        C0=transformed_coefficients,
        dt=0.002,
        steps=8,
        store_every=2,
        settings=settings,
    )
    first = run_self_consistent_block_dynamics_v214(
        provider(),
        provenance,
        initial_basis=initial_basis,
        C0=transformed_coefficients,
        dt=0.002,
        steps=3,
        store_every=1,
        settings=settings,
    )
    resumed = run_self_consistent_block_dynamics_v214(
        provider(),
        provenance,
        checkpoint=first["checkpoint"],
        steps=5,
        store_every=1,
        settings=settings,
    )

    assert max(
        np.linalg.norm(left.q - right.q)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) < 1.0e-13
    assert max(
        np.linalg.norm(left.p - right.p)
        for left, right in zip(full["final_basis"], resumed["final_basis"])
    ) < 1.0e-13
    assert _phase_aligned_metric_error(
        full["final_coefficients"], resumed["final_coefficients"], full["final_S"]
    ) < 2.0e-12
