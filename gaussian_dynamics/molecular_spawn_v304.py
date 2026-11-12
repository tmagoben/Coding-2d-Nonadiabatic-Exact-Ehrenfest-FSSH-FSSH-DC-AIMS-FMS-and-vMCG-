"""Replayable residual-oracle checkpoints from the complete molecular Hamiltonian.

Every generated candidate is evaluated, including rejections. These are local
residual-screening labels, not final topology decisions or certified ML records.
In particular, a full electronic spinor has no single target-state label.
"""
from dataclasses import asdict, replace
from functools import lru_cache
import json

from .correlated_basis_adaptation_v270 import (
    ControlledMultidimensionalBasisSettingsV270,
    generate_multidimensional_spawn_candidates_v270,
)
from .correlated_gaussian_tdvp_v270 import CorrelatedVariationalSettingsV270
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_geometry_v293 import straight_path_v293, wilson_path_v293
from .molecular_hamiltonian_v301 import build_molecular_tdvp_snapshot_v301
from .molecular_midpoint_v302 import _state_payload, _state_from_payload
from .pre_ml_contract_v294 import (
    V294_ORACLE_IMPLEMENTATION, evaluate_candidate_oracle_v294, stable_hash_v294,
)

SCHEMA_V304 = "gnd-molecular-spawn-checkpoint-v0.30.4-dev"


def _checkpoint_connection(bridge):
    """Cache repeated spinor columns locally; never reuse across checkpoints."""
    source = bridge.connection_field(transport_mode="nac_ode")

    @lru_cache(maxsize=None)
    def segment(start, end):
        return wilson_path_v293(source, straight_path_v293(start, end))

    return replace(source, segment_transport=lambda a, b: segment(tuple(a), tuple(b)).copy())


def _name(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _restore_settings(payload):
    values = dict(payload)
    inherited = dict(values.pop("inherited"))
    inherited["tdvp_settings"] = CorrelatedVariationalSettingsV270(**inherited["tdvp_settings"])
    return CurvedBasisSettingsV292(
        inherited=ControlledMultidimensionalBasisSettingsV270(**inherited), **values,
    ).validate()


def evaluate_molecular_spawn_checkpoint_v304(
    state, bridge, hamiltonian, *, family_id, trajectory_id, step,
    packet_ids=None, settings=CurvedBasisSettingsV292(), tangent_step=2e-6,
):
    """Rebuild H psi and TDVP, then record all deterministic candidate scores.

    The straight Wilson path and NAC-ODE transport match the v0.30.2 integrator.
    No state is changed. The package v0.29.4 split is retained only as historical
    metadata; no dataset split or training eligibility is certified here.
    """
    state.validate(require_normalized=False)
    settings.validate()
    family_id, trajectory_id = _name(family_id, "family_id"), _name(trajectory_id, "trajectory_id")
    if type(step) is not int or step < 0:
        raise ValueError("step must be a non-negative integer")
    ids = tuple(f"g{i:06d}" for i in range(state.ngaussian)) if packet_ids is None else tuple(packet_ids)
    if len(ids) != state.ngaussian or any(not isinstance(x, str) or not x.strip() for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("packet IDs must be non-empty, unique and complete")
    evaluation = build_molecular_tdvp_snapshot_v301(
        state, bridge, hamiltonian, relative_step=tangent_step,
    )
    connection = _checkpoint_connection(bridge)
    candidates = generate_multidimensional_spawn_candidates_v270(state, settings=settings.inherited)
    # JSON tuple encoding avoids delimiter collisions in family/trajectory IDs.
    group = json.dumps([family_id, trajectory_id], separators=(",", ":"), ensure_ascii=True)
    rows = []
    for serial, candidate in enumerate(candidates):
        record = evaluate_candidate_oracle_v294(
            state, connection, hamiltonian.quadrature, evaluation.snapshot, candidate,
            trajectory_group_id=group, packet_ids=ids, settings=settings,
        )
        rows.append({
            "serial": serial, "source_packet_id": ids[candidate.source_packet],
            "candidate": candidate.as_dict(), "record": record.as_dict(),
            "record_fingerprint": record.record_fingerprint,
        })
    receipt = {
        "schema": SCHEMA_V304, "oracle_implementation": V294_ORACLE_IMPLEMENTATION,
        "identity": {"family_id": family_id, "trajectory_id": trajectory_id, "step": step},
        "packet_ids": list(ids), "state": _state_payload(state),
        "hamiltonian_fingerprint": hamiltonian.fingerprint,
        "bridge_fingerprint": bridge.bridge_fingerprint,
        "settings": asdict(settings), "tangent_step": float(tangent_step),
        "candidate_count": len(rows), "candidates": rows,
        "candidate_space": "full_complex_electronic_spinor",
        "decision_stage": "residual_screening", "final_admission": None,
        "dataset_certified": False,
        "empty_candidate_reason": "no-nondegenerate-principal-axes" if not rows else None,
    }
    receipt["fingerprint"] = stable_hash_v294(receipt)
    return receipt


def serialize_molecular_spawn_checkpoint_v304(receipt):
    payload = dict(receipt)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != SCHEMA_V304 or stable_hash_v294(payload) != digest:
        raise ValueError("spawn checkpoint schema or fingerprint mismatch")
    return json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"


def replay_molecular_spawn_checkpoint_v304(serialized, bridge, hamiltonian):
    """Recompute every label; hashes alone do not establish scientific correctness.

    Exact replay requires the same numerical runtime, provider and Hamiltonian.
    Gauge comparisons require tolerances and are a separate physical audit.
    """
    receipt = json.loads(serialized)
    serialize_molecular_spawn_checkpoint_v304(receipt)
    if receipt["hamiltonian_fingerprint"] != hamiltonian.fingerprint:
        raise ValueError("replay Hamiltonian fingerprint mismatch")
    rebuilt = evaluate_molecular_spawn_checkpoint_v304(
        _state_from_payload(receipt["state"]), bridge, hamiltonian,
        **receipt["identity"], packet_ids=receipt["packet_ids"],
        settings=_restore_settings(receipt["settings"]), tangent_step=receipt["tangent_step"],
    )
    if rebuilt != receipt:
        raise ValueError("spawn replay does not reproduce stored checkpoint")
    return rebuilt
