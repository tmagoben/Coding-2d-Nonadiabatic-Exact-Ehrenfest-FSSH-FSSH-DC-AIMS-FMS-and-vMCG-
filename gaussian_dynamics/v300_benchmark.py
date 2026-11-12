"""Deterministic evidence campaign for the v0.30.0 molecular-frame bridge."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220
from .complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .molecular_frame_bridge_v300 import (
    MolecularFrameBridgeSettingsV300,
    MolecularFrameBridgeV300,
    molecular_spawn_context_v300,
    serialize_molecular_spawn_contexts_v300,
)
from .pre_ml_contract_v294 import (
    SpawnOracleRecordV294,
    V294_DATASET_SCHEMA,
    V294_ORACLE_SCHEMA,
    canonical_state_payload_v294,
    stable_hash_v294,
)


V300_EVIDENCE_SCHEMA = "gnd-molecular-frame-bridge-evidence-v0.30.0-dev"


def _fixture():
    base = AnalyticSingletTripletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        U0=random_unitary_v21(4, seed=3000),
        phase_gradient=np.asarray([[0.18], [-0.11], [0.07], [0.23]]),
        phase_offset=np.asarray([0.12, -0.08, 0.05, -0.03]),
    )
    moving = GaugeTransformedOperatorProviderV21(base, gauge)
    bridge = MolecularFrameBridgeV300(
        moving,
        np.asarray([0.0]),
        settings=MolecularFrameBridgeSettingsV300(
            connection_overlap_tolerance=3.0e-9,
            maximum_audit_segment_length=0.6,
        ),
    )
    return base, gauge, moving, bridge


def _context_fixture(bridge):
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
        V294_ORACLE_SCHEMA,
        V294_DATASET_SCHEMA,
        "trajectory-v300-fixture",
        "validation",
        state_fingerprint,
        stable_hash_v294({"candidate": "fixture"}),
        stable_hash_v294({"oracle": "fixture"}),
        {"feature_schema": "gnd-spawn-features-invariant-v1"},
        0.02,
        0.1,
        2,
        3.0,
        (1.0, 0.5),
        True,
        (),
    ).validate()
    return molecular_spawn_context_v300(
        record,
        state,
        bridge,
        packet_ids=packet_ids,
        path_family=path_family,
    )


def run_v300_benchmark(output_path=None):
    base, gauge, moving, bridge = _fixture()
    q = np.asarray([0.17])
    frame = gauge.matrix(q)
    fixed = base.evaluate(q)
    transformed = bridge.operator_point(q)
    h_covariance = float(
        np.linalg.norm(transformed.H - frame.conj().T @ fixed.H @ frame, ord="fro")
    )
    k_covariance = float(
        np.linalg.norm(
            transformed.dH_dq[0] - frame.conj().T @ fixed.dH_dq[0] @ frame,
            ord="fro",
        )
    )
    connection_covariance = float(
        np.linalg.norm(
            transformed.connection_q[0]
            - frame.conj().T @ gauge.derivatives(q)[0],
            ord="fro",
        )
    )

    receipts = bridge.audit_path(
        (np.asarray([-0.18]), np.asarray([0.02]), np.asarray([0.24]))
    )
    maximum_transport_residual = float(
        max(receipt.connection_overlap_residual for receipt in receipts)
    )
    maximum_overlap_unitarity = float(
        max(receipt.overlap_transport_unitarity_residual for receipt in receipts)
    )
    maximum_connection_unitarity = float(
        max(receipt.connection_transport_unitarity_residual for receipt in receipts)
    )

    context = _context_fixture(bridge)
    dataset_once = serialize_molecular_spawn_contexts_v300([context])
    dataset_twice = serialize_molecular_spawn_contexts_v300([context])

    class ExpandedOverlapProvider:
        provenance = base.provenance

        def evaluate_snapshot(self, point):
            return moving.evaluate_snapshot(point)

        def snapshot_overlap(self, left, right):
            return 1.02 * moving.snapshot_overlap(left, right)

    raw_overlap_rejected = False
    try:
        bad = MolecularFrameBridgeV300(
            ExpandedOverlapProvider(), np.asarray([0.0])
        )
        bad.overlap_transport(np.asarray([-0.1]), np.asarray([0.1]))
    except ValueError:
        raw_overlap_rejected = True

    strict_identity_rejected = False
    try:
        MolecularFrameBridgeV300(
            moving,
            np.asarray([0.0]),
            settings=MolecularFrameBridgeSettingsV300(
                require_numerical_identity=True
            ),
        )
    except ValueError:
        strict_identity_rejected = True

    metrics = {
        "hamiltonian_gauge_covariance_residual": h_covariance,
        "operator_gradient_gauge_covariance_residual": k_covariance,
        "connection_gauge_covariance_residual": connection_covariance,
        "maximum_connection_overlap_transport_residual": maximum_transport_residual,
        "maximum_overlap_transport_unitarity_residual": maximum_overlap_unitarity,
        "maximum_connection_transport_unitarity_residual": maximum_connection_unitarity,
        "molecular_context_fingerprint": context.fingerprint,
        "provider_fingerprint": bridge.provider_fingerprint,
        "bridge_fingerprint": bridge.bridge_fingerprint,
        "dataset_sha256": stable_hash_v294(dataset_once),
    }
    gates = {
        "total_soc_hamiltonian_is_gauge_covariant": h_covariance < 5.0e-12,
        "physical_operator_gradient_is_gauge_covariant": k_covariance < 5.0e-12,
        "derivative_connection_matches_gauge_rule": connection_covariance < 5.0e-12,
        "nac_wilson_and_overlap_polar_transport_agree": maximum_transport_residual < 3.0e-9,
        "overlap_polar_transport_is_unitary": maximum_overlap_unitarity < 5.0e-12,
        "connection_wilson_transport_is_unitary": maximum_connection_unitarity < 5.0e-10,
        "raw_overlap_spectral_expansion_fails_closed": raw_overlap_rejected,
        "production_numerical_identity_gate_fails_closed": strict_identity_rejected,
        "molecular_context_preserves_group_split": context.split == "validation",
        "molecular_context_retains_soc_provenance": context.soc_enabled,
        "molecular_context_serialization_is_deterministic": dataset_once == dataset_twice,
        "molecular_context_excludes_raw_frame_arrays": all(
            token not in context.as_dict()
            for token in ("H", "connection_q", "coefficients")
        ),
    }
    payload = {
        "schema": V300_EVIDENCE_SCHEMA,
        "classification": "development checkpoint; no molecular trajectory claim",
        "metrics": metrics,
        "gates": {name: bool(value) for name, value in gates.items()},
        "passed": bool(all(gates.values())),
    }
    payload["fingerprint"] = stable_hash_v294(payload)
    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    return payload


if __name__ == "__main__":
    evidence = run_v300_benchmark(
        Path(__file__).resolve().parents[1]
        / "results"
        / "v0300_molecular_frame_bridge_evidence.json"
    )
    print(json.dumps(evidence, indent=2, sort_keys=True, allow_nan=False))
