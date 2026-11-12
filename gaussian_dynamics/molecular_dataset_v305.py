"""Explicit full-spinor admission adapter; replay verification is not certification.

The frozen, state-targeted v0.29.4 schema cannot represent these events. This
version keeps geometry features, oracle targets and provenance in separate
fields, with no fabricated target state or training eligibility.
"""
from dataclasses import dataclass
import json
import re

import numpy as np

from .molecular_admission_v305 import replay_molecular_admission_v305
from .pre_ml_contract_v294 import stable_hash_v294

DATASET_SCHEMA_V305 = "gnd-full-spinor-admission-dataset-v0.30.5-dev"
SPLITS_V305 = ("train", "validation", "test", "audit")
# These are scalar geometry/row-norm features. Parameter-metric diagnostics and
# residual scores are deliberately excluded from the predictor's inputs.
FEATURE_KEYS_V305 = (
    "ngaussian", "ndim", "nstate", "dq_width_norm2", "dp_width_inverse_norm2",
    "dq_euclidean_norm", "dp_euclidean_norm", "width_eigenvalues",
    "chirp_eigenvalues", "width_logdet", "chirp_frobenius_norm", "displacement_kind",
)


@dataclass(frozen=True)
class MolecularAdmissionSourceV305:
    serialized: str
    bridge: object
    hamiltonian: object
    variant_id: str = "canonical"


def trajectory_split_v305(family_id, trajectory_id, *, salt="v0.30.5", fractions=(.70, .15, .10, .05)):
    """Tuple encoding prevents delimiter collisions; every variant shares a split."""
    if any(not isinstance(x, str) or not x.strip() for x in (family_id, trajectory_id, salt)):
        raise ValueError("split identities and salt must be nonempty strings")
    fractions = tuple(float(x) for x in fractions)
    if len(fractions) != 4 or not all(np.isfinite(x) and x >= 0 for x in fractions) or abs(sum(fractions) - 1) > 1e-12:
        raise ValueError("four finite nonnegative split fractions must sum to one")
    group = [family_id, trajectory_id]
    integer = int(stable_hash_v294({"salt": salt, "group": group})[:16], 16)
    # Compare integers with cuts, avoiding a rounded u == 1 at the top bin.
    cumulative = 0.0
    for name, fraction in zip(SPLITS_V305[:-1], fractions[:-1]):
        cumulative += fraction
        if integer < int(cumulative * 2**64):
            return name
    return SPLITS_V305[-1]


def materialize_admission_dataset_v305(sources, *, generation_commit, salt="v0.30.5", fractions=(.70, .15, .10, .05)):
    """Replay every checkpoint before exporting every candidate, including blocks.

    Final admission is a context-dependent selection target. A candidate marked
    not_selected/blocked has not been found physically inadmissible. The manifest
    stays uncertified until a separate physical symmetry campaign is implemented.
    """
    if not isinstance(generation_commit, str) or re.fullmatch(r"[0-9a-f]{40}", generation_commit) is None:
        raise ValueError("generation_commit must be a full lowercase Git SHA")
    sources = tuple(sources)
    if not sources:
        raise ValueError("dataset requires at least one replayable checkpoint")
    fractions = tuple(fractions)
    # Validate split controls before performing any expensive physical replay.
    trajectory_split_v305("validation", "validation", salt=salt, fractions=fractions)
    rows, checkpoints, seen = [], [], set()
    for source in sources:
        variant = source.variant_id
        if not isinstance(variant, str) or not variant.strip() or variant != variant.strip():
            raise ValueError("variant ID must be a nonempty trimmed string")
        receipt = replay_molecular_admission_v305(source.serialized, source.bridge, source.hamiltonian).receipt
        screening = receipt["screening"]
        identity = screening["identity"]
        key = (identity["family_id"], identity["trajectory_id"], identity["step"], variant)
        if key in seen:
            raise ValueError("duplicate checkpoint identity and variant")
        seen.add(key)
        split = trajectory_split_v305(identity["family_id"], identity["trajectory_id"], salt=salt, fractions=fractions)
        checkpoints.append({"identity": identity, "variant_id": variant, "split": split,
                            "fingerprint": receipt["fingerprint"], "candidate_count": screening["candidate_count"]})
        for candidate, decision in zip(screening["candidates"], receipt["candidate_decisions"]):
            native = candidate["record"]
            features = {name: native["features"][name] for name in FEATURE_KEYS_V305}
            features["source_coefficient_norm_squared"] = native["features"]["source_population"]
            row = {
                "schema": DATASET_SCHEMA_V305, "identity": {**identity, "variant_id": variant, "candidate_serial": candidate["serial"], "source_packet_id": candidate["source_packet_id"]},
                "split": split, "candidate_space": receipt["candidate_space"], "target_state": None,
                "features": features,
                "targets": {"predicted_residual_reduction": native["predicted_residual_reduction"],
                            "relative_residual_reduction": native["relative_residual_reduction"],
                            "residual_accepted": decision["residual_accepted"],
                            "residual_rejection_reasons": native["rejection_reasons"],
                            "final_admitted": decision["final_admitted"], "outcome": decision["outcome"],
                            "projection_attempted": decision["attempt"] is not None},
                "provenance": {"admission_fingerprint": receipt["fingerprint"],
                               "screening_fingerprint": screening["fingerprint"],
                               "oracle_record_fingerprint": candidate["record_fingerprint"],
                               "oracle_implementation": screening["oracle_implementation"],
                               "generation_commit": generation_commit},
            }
            row["fingerprint"] = stable_hash_v294(row)
            rows.append(row)
    rows.sort(key=lambda r: (r["identity"]["family_id"], r["identity"]["trajectory_id"], r["identity"]["step"], r["identity"]["variant_id"], r["identity"]["candidate_serial"]))
    checkpoints.sort(key=lambda c: (c["identity"]["family_id"], c["identity"]["trajectory_id"], c["identity"]["step"], c["variant_id"]))
    payload = {
        "schema": DATASET_SCHEMA_V305, "generation_commit": generation_commit,
        "split_policy": {"salt": salt, "fractions": list(fractions), "group": "family_id + trajectory_id; all steps and variants"},
        "checkpoints": checkpoints, "rows": rows,
        "checkpoint_count": len(checkpoints), "record_count": len(rows),
        "split_counts": {name: sum(row["split"] == name for row in rows) for name in SPLITS_V305},
        "exact_replay_verified": True, "physical_invariance_certified": False,
        "dataset_certified": False, "training_eligible": False,
    }
    payload["fingerprint"] = stable_hash_v294(payload)
    return payload


def admission_dataset_jsonl_v305(dataset):
    """Export sealed rows; consumers must retain the accompanying manifest."""
    payload = dict(dataset)
    digest = payload.pop("fingerprint")
    if payload.get("schema") != DATASET_SCHEMA_V305 or stable_hash_v294(payload) != digest:
        raise ValueError("dataset schema or fingerprint mismatch")
    return "".join(json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n" for row in dataset["rows"])
