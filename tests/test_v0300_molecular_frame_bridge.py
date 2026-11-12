import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticSingletTripletSOCProviderV220,
)
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
)
from gaussian_dynamics.curved_geometry_v293 import (
    straight_path_v293,
    wilson_path_v293,
)
from gaussian_dynamics.molecular_frame_bridge_v300 import (
    MolecularFrameBridgeSettingsV300,
    MolecularFrameBridgeV300,
    molecular_spawn_context_v300,
    serialize_molecular_spawn_contexts_v300,
)
from gaussian_dynamics.pre_ml_contract_v294 import (
    SpawnOracleRecordV294,
    V294_DATASET_SCHEMA,
    V294_ORACLE_SCHEMA,
    canonical_state_payload_v294,
    stable_hash_v294,
)
from gaussian_dynamics.v300_benchmark import run_v300_benchmark


def moving_provider_fixture():
    base = AnalyticSingletTripletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        U0=random_unitary_v21(4, seed=3000),
        phase_gradient=np.asarray([[0.18], [-0.11], [0.07], [0.23]]),
        phase_offset=np.asarray([0.12, -0.08, 0.05, -0.03]),
    )
    return base, gauge, GaugeTransformedOperatorProviderV21(base, gauge)


def bridge_fixture():
    _, _, provider = moving_provider_fixture()
    settings = MolecularFrameBridgeSettingsV300(
        connection_overlap_tolerance=3.0e-9,
        maximum_audit_segment_length=0.6,
    )
    return MolecularFrameBridgeV300(provider, np.asarray([0.0]), settings=settings)


def test_moving_molecular_provider_exposes_total_soc_operator_and_nac_connection():
    base, gauge, provider = moving_provider_fixture()
    bridge = MolecularFrameBridgeV300(provider, np.asarray([0.0]))
    q = np.asarray([0.17])
    fixed = base.evaluate(q)
    moving = bridge.operator_point(q)
    frame = gauge.matrix(q)
    assert bridge.provenance.soc_enabled
    assert np.linalg.norm(moving.H - frame.conj().T @ fixed.H @ frame) < 2.0e-13
    assert np.linalg.norm(
        moving.dH_dq[0] - frame.conj().T @ fixed.dH_dq[0] @ frame
    ) < 2.0e-13
    assert np.linalg.norm(
        moving.connection_q[0] - frame.conj().T @ gauge.derivatives(q)[0]
    ) < 2.0e-13
    assert np.linalg.norm(moving.H - moving.H.conj().T) < 2.0e-13
    assert np.linalg.norm(
        moving.connection_q[0] + moving.connection_q[0].conj().T
    ) < 2.0e-13


def test_nac_wilson_transport_matches_certified_overlap_polar_transport():
    bridge = bridge_fixture()
    start = np.asarray([-0.18])
    end = np.asarray([0.24])
    receipt = bridge.audit_segment(start, end)
    assert receipt.passed
    assert receipt.connection_overlap_residual < 3.0e-9
    assert receipt.overlap_reciprocity_residual < 2.0e-12
    assert receipt.transport_reciprocity_residual < 2.0e-12
    assert receipt.overlap_transport_unitarity_residual < 2.0e-12
    assert receipt.connection_transport_unitarity_residual < 2.0e-10


def test_overlap_polar_mode_is_explicit_and_matches_direct_certified_transport():
    bridge = bridge_fixture()
    start = np.asarray([-0.12])
    end = np.asarray([0.19])
    field = bridge.connection_field(transport_mode="overlap_polar")
    transported = wilson_path_v293(field, straight_path_v293(start, end))
    expected = bridge.overlap_transport(start, end)
    assert np.linalg.norm(transported - expected) < 2.0e-13
    with pytest.raises(ValueError, match="transport_mode"):
        bridge.connection_field(transport_mode="raw-overlap")


def test_electronic_only_action_is_gauge_covariant():
    base, gauge, provider = moving_provider_fixture()
    bridge = MolecularFrameBridgeV300(provider, np.asarray([0.0]))
    points = np.asarray([[-0.2], [0.0], [0.27]])
    rng = np.random.default_rng(3001)
    fixed_vectors = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
    moving_vectors = np.asarray(
        [gauge.matrix(q).conj().T @ vector for q, vector in zip(points, fixed_vectors)]
    )
    observed = bridge.apply_electronic_hamiltonian(points, moving_vectors)
    expected = np.asarray(
        [
            gauge.matrix(q).conj().T @ (base.evaluate(q).H @ vector)
            for q, vector in zip(points, fixed_vectors)
        ]
    )
    assert np.max(np.abs(observed - expected)) < 3.0e-13


def test_raw_spectral_expansion_is_rejected_before_transport():
    _, _, provider = moving_provider_fixture()

    class ExpandedOverlapProvider:
        provenance = provider.base_provider.provenance

        def evaluate_snapshot(self, q):
            return provider.evaluate_snapshot(q)

        def snapshot_overlap(self, left, right):
            return 1.02 * provider.snapshot_overlap(left, right)

    bridge = MolecularFrameBridgeV300(ExpandedOverlapProvider(), np.asarray([0.0]))
    with pytest.raises(ValueError, match="physically inconsistent"):
        bridge.overlap_transport(np.asarray([-0.1]), np.asarray([0.1]))


def test_nonantihermitian_derivative_connection_fails_closed():
    _, _, provider = moving_provider_fixture()

    class BrokenConnectionProvider:
        provenance = provider.base_provider.provenance

        def evaluate_snapshot(self, q):
            snapshot = provider.evaluate_snapshot(q)
            snapshot.point.connection_q[0, 0, 0] += 0.1
            return snapshot

        def snapshot_overlap(self, left, right):
            return provider.snapshot_overlap(left, right)

    with pytest.raises(ValueError, match="anti-Hermiticity"):
        MolecularFrameBridgeV300(BrokenConnectionProvider(), np.asarray([0.0]))


def test_production_identity_gate_rejects_provenance_only_provider():
    _, _, provider = moving_provider_fixture()
    settings = MolecularFrameBridgeSettingsV300(require_numerical_identity=True)
    with pytest.raises(ValueError, match="numerical identity"):
        MolecularFrameBridgeV300(provider, np.asarray([0.0]), settings=settings)


def test_molecular_spawn_context_preserves_group_split_and_excludes_raw_frames():
    bridge = bridge_fixture()
    state = CorrelatedGaussianSpinorStateV270(
        q=np.asarray([[0.13]]),
        p=np.asarray([[0.07]]),
        width_matrices=np.asarray([[[1.1]]]),
        chirp_matrices=np.asarray([[[0.03]]]),
        coefficients=np.asarray([[0.72 + 0.03j, 0.31 - 0.08j, 0.2j, 0.11]]),
    ).normalized()
    packet_ids = ("packet-A",)
    path_family = "straight_path_v293"
    state_fingerprint = stable_hash_v294(
        canonical_state_payload_v294(
            state, packet_ids=packet_ids, path_family=path_family
        )
    )
    record = SpawnOracleRecordV294(
        oracle_schema=V294_ORACLE_SCHEMA,
        dataset_schema=V294_DATASET_SCHEMA,
        trajectory_group_id="trajectory-v300-fixture",
        split="validation",
        state_fingerprint=state_fingerprint,
        candidate_fingerprint=stable_hash_v294({"candidate": "fixture"}),
        oracle_config_fingerprint=stable_hash_v294({"oracle": "fixture"}),
        features={"feature_schema": "gnd-spawn-features-invariant-v1"},
        predicted_residual_reduction=0.02,
        relative_residual_reduction=0.1,
        schur_rank=2,
        schur_condition_number=3.0,
        schur_singular_values=(1.0, 0.5),
        accepted_by_residual_oracle=True,
        rejection_reasons=(),
    ).validate()
    context = molecular_spawn_context_v300(
        record,
        state,
        bridge,
        packet_ids=packet_ids,
        path_family=path_family,
    )
    assert context.split == record.split
    assert context.trajectory_group_id == record.trajectory_group_id
    assert context.state_fingerprint == record.state_fingerprint
    assert context.soc_enabled
    payload = context.as_dict()
    assert "connection_q" not in payload
    assert "coefficients" not in payload
    assert "H" not in payload
    first = serialize_molecular_spawn_contexts_v300([context])
    second = serialize_molecular_spawn_contexts_v300([context])
    assert first == second
    assert context.fingerprint == context.fingerprint


def test_molecular_spawn_context_rejects_state_or_path_mismatch():
    bridge = bridge_fixture()
    state = CorrelatedGaussianSpinorStateV270(
        q=np.asarray([[0.0]]),
        p=np.asarray([[0.0]]),
        width_matrices=np.asarray([[[1.0]]]),
        chirp_matrices=np.asarray([[[0.0]]]),
        coefficients=np.asarray([[1.0, 0.0, 0.0, 0.0]], dtype=complex),
    ).validate(require_normalized=False)
    record = SpawnOracleRecordV294(
        V294_ORACLE_SCHEMA,
        V294_DATASET_SCHEMA,
        "group",
        "train",
        stable_hash_v294({"different": "state"}),
        stable_hash_v294({"candidate": 1}),
        stable_hash_v294({"settings": 1}),
        {},
        0.0,
        0.0,
        0,
        1.0,
        (),
        False,
        ("no-new-curved-tangent-direction",),
    ).validate()
    with pytest.raises(ValueError, match="does not match"):
        molecular_spawn_context_v300(record, state, bridge)


def test_v300_evidence_campaign_passes_all_native_gates():
    evidence = run_v300_benchmark()
    assert evidence["schema"] == "gnd-molecular-frame-bridge-evidence-v0.30.0-dev"
    assert evidence["passed"] is True
    assert len(evidence["gates"]) == 12
    assert all(type(value) is bool for value in evidence["gates"].values())
