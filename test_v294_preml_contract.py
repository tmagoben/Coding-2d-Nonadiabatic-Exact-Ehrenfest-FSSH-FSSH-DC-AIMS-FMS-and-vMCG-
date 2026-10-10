import numpy as np
import pytest

from v294_preml_contract import (
    DatasetManifest,
    EventIdentity,
    InvariantFeatures,
    OracleDecision,
    SCHEMA_VERSION,
    SpawningRecord,
    assert_disjoint_manifests,
    canonical_json,
    evaluate_predictions,
    fingerprint,
    hermitian_operator_invariants,
)


def make_record(*, trajectory="traj-a", family="model-a", serial=0, action="spawn", gain=0.3):
    identity = EventIdentity(
        trajectory_id=trajectory,
        family_id=family,
        event_id=f"event-{serial}",
        packet_id="p0",
        parent_packet_id=None,
        step=10 + serial,
        time_au=0.5 + 0.1 * serial,
    )
    features = InvariantFeatures(
        populations=(0.7, 0.3),
        energy_gaps=(0.015,),
        force_norms=(0.02, 0.03),
        nac_frobenius=0.12,
        soc_singular_values=(0.001, 0.003),
        berry_curvature_frobenius=0.05,
        gram_min_eigenvalue=0.2,
        gram_condition_number=8.0,
        max_packet_overlap_abs=0.4,
        nuclear_speed=0.08,
        kinetic_energy=0.02,
        residual_norm=0.9,
        novelty=0.6,
        projection_loss_estimate=1e-7,
        dt_au=0.01,
    )
    decision = OracleDecision(
        action=action,
        target_state=1 if action == "spawn" else None,
        residual_gain=gain,
        normalized_residual_gain=abs(gain),
        threshold=0.1,
        reason_code="residual_gain_pass" if action == "spawn" else "below_threshold",
        oracle_version="v0.29.3-oracle",
        admitted=True if action == "spawn" else False,
    )
    return SpawningRecord(
        identity=identity,
        features=features,
        decision=decision,
        candidate_kind="principal_axis_displacement",
        candidate_serial=serial,
    )


def test_schema_version_is_frozen():
    assert SCHEMA_VERSION == "0.29.4"


def test_canonical_json_and_fingerprint_ignore_mapping_order():
    left = {"b": 2, "a": {"y": 4, "x": 3}}
    right = {"a": {"x": 3, "y": 4}, "b": 2}
    assert canonical_json(left) == canonical_json(right)
    assert fingerprint(left) == fingerprint(right)


def test_record_hash_is_deterministic_and_physics_sensitive():
    a = make_record(gain=0.3)
    b = make_record(gain=0.3)
    c = make_record(gain=0.31)
    assert a.record_hash == b.record_hash
    assert a.record_hash != c.record_hash


def test_nonfinite_feature_is_rejected():
    with pytest.raises(ValueError, match="finite"):
        InvariantFeatures(populations=(0.5, 0.5), residual_norm=float("nan"))


def test_populations_must_be_normalized():
    with pytest.raises(ValueError, match="sum to 1"):
        InvariantFeatures(populations=(0.8, 0.3))


def test_spawn_and_no_spawn_target_state_contract():
    with pytest.raises(ValueError, match="target_state"):
        OracleDecision(
            action="spawn",
            target_state=None,
            residual_gain=0.1,
            normalized_residual_gain=0.1,
            threshold=0.05,
            reason_code="pass",
            oracle_version="oracle",
        )
    with pytest.raises(ValueError, match="must not carry"):
        OracleDecision(
            action="no_spawn",
            target_state=1,
            residual_gain=0.0,
            normalized_residual_gain=0.0,
            threshold=0.05,
            reason_code="reject",
            oracle_version="oracle",
        )


def test_native_boolean_evidence_boundary_is_enforced():
    with pytest.raises(TypeError, match="native bool"):
        OracleDecision(
            action="spawn",
            target_state=1,
            residual_gain=0.2,
            normalized_residual_gain=0.2,
            threshold=0.1,
            reason_code="pass",
            oracle_version="oracle",
            admitted=np.bool_(True),
        )


def test_dataset_hash_is_independent_of_record_order():
    records = [
        make_record(trajectory="t1", serial=1),
        make_record(trajectory="t2", serial=2, action="no_spawn", gain=0.02),
    ]
    a = DatasetManifest.from_records("train", records, "abc123")
    b = DatasetManifest.from_records("train", reversed(records), "abc123")
    assert a.dataset_hash == b.dataset_hash
    assert a.record_hashes == b.record_hashes


def test_leakage_guard_rejects_same_family_trajectory_across_splits():
    train = DatasetManifest.from_records("train", [make_record(trajectory="same")], "abc123")
    test = DatasetManifest.from_records(
        "test", [make_record(trajectory="same", serial=2, action="no_spawn")], "abc123"
    )
    with pytest.raises(ValueError, match="leakage group"):
        assert_disjoint_manifests(train, test)


def test_leakage_guard_allows_different_trajectories():
    train = DatasetManifest.from_records("train", [make_record(trajectory="t1")], "abc123")
    test = DatasetManifest.from_records("test", [make_record(trajectory="t2", serial=2)], "abc123")
    assert_disjoint_manifests(train, test)


def test_hermitian_invariants_survive_unitary_similarity():
    h = np.array([[0.1, 0.02 + 0.03j], [0.02 - 0.03j, 0.4]], dtype=complex)
    raw = np.array([[1 + 1j, 2 - 0.3j], [0.7 + 0.2j, -1 + 0.5j]], dtype=complex)
    q, r = np.linalg.qr(raw)
    phases = np.diag(np.exp(-1j * np.angle(np.diag(r))))
    u = q @ phases
    transformed = u.conj().T @ h @ u
    a = hermitian_operator_invariants(h)
    b = hermitian_operator_invariants(transformed)
    assert np.isclose(a["trace"], b["trace"], atol=1e-12)
    assert np.isclose(a["frobenius_norm"], b["frobenius_norm"], atol=1e-12)
    assert np.allclose(a["eigenvalues"], b["eigenvalues"], atol=1e-12)
    assert np.allclose(a["adjacent_gaps"], b["adjacent_gaps"], atol=1e-12)


def test_nonhermitian_operator_is_rejected():
    with pytest.raises(ValueError, match="Hermitian"):
        hermitian_operator_invariants(np.array([[0.0, 1.0], [0.0, 0.0]], dtype=complex))


def test_evaluation_exposes_false_negative_rate_and_score_error():
    truth = [make_record(serial=0, action="spawn", gain=0.3), make_record(serial=1, action="no_spawn", gain=0.02)]
    result = evaluate_predictions(truth, probabilities=[0.4, 0.2], predicted_scores=[0.25, 0.03])
    assert result["false_negative_rate"] == 1.0
    assert result["recall"] == 0.0
    assert np.isclose(result["mean_score_error"], 0.03)


def test_manifest_payload_is_json_serializable_and_self_consistent():
    manifest = DatasetManifest.from_records("audit", [make_record()], "deadbeef")
    payload = manifest.to_payload()
    assert payload["record_count"] == 1
    assert payload["schema_version"] == SCHEMA_VERSION
    canonical_json(payload)
