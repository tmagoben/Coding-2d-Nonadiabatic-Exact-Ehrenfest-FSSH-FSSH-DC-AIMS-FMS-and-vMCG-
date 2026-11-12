"""Audit complete-complex spawning scores on actual finite-grid SOC residuals."""
from dataclasses import replace
import platform

import numpy as np
import scipy

from .correlated_basis_adaptation_v270 import MultidimensionalSpawnCandidateV270
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import candidate_coefficient_block_v293
from .molecular_hamiltonian_v301 import build_molecular_tdvp_snapshot_v301
from .molecular_midpoint_v302 import molecular_midpoint_step_v302
from .molecular_spawn_v304 import (
    evaluate_molecular_spawn_checkpoint_v304, replay_molecular_spawn_checkpoint_v304,
    serialize_molecular_spawn_checkpoint_v304, _checkpoint_connection,
)
from .pre_ml_contract_v294 import stable_hash_v294
from .v303_benchmark import soc_packet_fixture_v303

SCHEMA_V304_EVIDENCE = "gnd-molecular-spawn-audit-v0.30.4-dev"


def _direct_gain(j, k, residual, *, complete=True):
    # Independently solve both least-squares systems, without the production
    # projected-block SVD or its stored predicted gain.
    a = np.column_stack([j, k, 1j*k] if complete else [j, k])
    jr, ar = np.vstack([j.real, j.imag]), np.vstack([a.real, a.imag])
    rr = np.r_[residual.real, residual.imag]
    before = rr - jr @ np.linalg.lstsq(jr, rr, rcond=1e-10)[0]
    after = rr - ar @ np.linalg.lstsq(ar, rr, rcond=1e-10)[0]
    return float(before @ before - after @ after)


def _references(state, bridge, grid, checkpoint):
    evaluation = build_molecular_tdvp_snapshot_v301(state, bridge, grid)
    snap = evaluation.snapshot
    connection = _checkpoint_connection(bridge)
    references = []
    for row in checkpoint["candidates"]:
        candidate = MultidimensionalSpawnCandidateV270(**row["candidate"])
        k = candidate_coefficient_block_v293(candidate, state.nstate, connection, grid.quadrature)
        references.append({
            "complete_lstsq_gain": _direct_gain(snap.jacobian, k, snap.residual_weighted),
            "legacy_real_only_gain": _direct_gain(snap.jacobian, k, snap.residual_weighted, complete=False),
        })
    electronic = np.einsum("nij,nj->ni", grid.electronic_blocks, evaluation.psi)
    kinetic_norm = float(np.linalg.norm(evaluation.hpsi-electronic))
    return references, kinetic_norm


def evaluate_spawn_gates_v304(evidence):
    """Recompute gates from measurements, never trusting cached pass flags."""
    points = evidence["snapshots"]
    rows = [row for point in points for row in point["comparisons"]]
    def near(a, b):
        return bool(np.isfinite(a) and np.isfinite(b) and abs(a-b) <= 3e-10 + 3e-6*max(abs(a), abs(b)))
    rejected = evidence["rejected_checkpoint"]["candidates"]
    def records_match(point):
        for frame in ("fixed", "moving"):
            checkpoint = point[frame]
            if len(checkpoint["candidates"]) != len(point["comparisons"]):
                return False
            payload = dict(checkpoint)
            digest = payload.pop("fingerprint")
            if stable_hash_v294(payload) != digest:
                return False
            for item, comparison in zip(checkpoint["candidates"], point["comparisons"]):
                record = item["record"]
                fields = {"gain": "predicted_residual_reduction", "relative_gain": "relative_residual_reduction",
                          "accepted": "accepted_by_residual_oracle", "reasons": "rejection_reasons", "rank": "schur_rank"}
                if item["serial"] != comparison["serial"] or any(comparison[frame+"_"+key] != record[value] for key, value in fields.items()):
                    return False
        return all(a["candidate"] == b["candidate"] for a, b in zip(point["fixed"]["candidates"], point["moving"]["candidates"]))
    gates = {
        "initial_and_propagated_snapshots": len(points) == 2 and points[0]["time_au"] == 0.0 and points[1]["time_au"] == 0.05,
        "complete_candidate_coverage": bool(rows) and all(
            point["fixed"]["candidate_count"] == point["moving"]["candidate_count"] == len(point["comparisons"]) == 4
            and len(point["fixed"]["candidates"]) == len(point["moving"]["candidates"]) == 4 for point in points),
        "full_kinetic_action": all(np.isfinite(point["kinetic_action_norm"]) and point["kinetic_action_norm"] > 1e-3 for point in points),
        "comparison_records_match_checkpoints": bool(points) and all(records_match(point) for point in points),
        "bounded_residual_gains": bool(rows) and all(
            np.isfinite(row[frame+"_gain"]) and row[frame+"_gain"] >= 0
            and np.isfinite(row[frame+"_relative_gain"]) and 0 <= row[frame+"_relative_gain"] <= 1+1e-10
            for row in rows for frame in ("fixed", "moving")),
        "independent_complex_lstsq": bool(rows) and all(near(row[frame+"_gain"], row[frame+"_reference_gain"]) for row in rows for frame in ("fixed", "moving")),
        "physical_gauge_scores": bool(rows) and all(near(row["fixed_gain"], row["moving_gain"]) and near(row["fixed_relative_gain"], row["moving_relative_gain"]) for row in rows),
        "physical_gauge_decisions": bool(rows) and all(row["fixed_accepted"] is row["moving_accepted"] and row["fixed_reasons"] == row["moving_reasons"] for row in rows),
        "physical_gauge_ranks": bool(rows) and all(row["fixed_rank"] == row["moving_rank"] and row["fixed_rank"] > 0 for row in rows),
        "imaginary_directions_change_scores": bool(rows) and max(abs(row["fixed_gain"]-row["fixed_legacy_gain"]) for row in rows) > 1e-7,
        "same_runtime_exact_replay": evidence["same_runtime_exact_replay"] is True,
        "residual_rejections_recorded": len(rejected) == 4 and all(
            row["record"]["accepted_by_residual_oracle"] is False
            and "insufficient-predicted-residual-reduction" in row["record"]["rejection_reasons"] for row in rejected),
        "residual_labels_are_not_final_admission": all(
            point[frame]["decision_stage"] == "residual_screening"
            and point[frame]["final_admission"] is None
            and point[frame]["dataset_certified"] is False
            for point in points for frame in ("fixed", "moving")),
    }
    return {name: bool(value) for name, value in gates.items()}


def run_v304_spawn_campaign(*, verbose=False):
    state, bridge, grid, gauge, _ = soc_packet_fixture_v303()
    _, moving_bridge, moving_grid, _, _ = soc_packet_fixture_v303(moving=True)
    step = molecular_midpoint_step_v302(state, bridge, grid, 0.05)
    points = []
    for index, current in enumerate((state, step.after)):
        if verbose:
            print(f"Scoring all candidates at t={current.time_au} in both electronic gauges", flush=True)
        # Compare the SAME physical state. Independently propagated gauges can
        # differ by time-discretization error and would confound this local audit.
        transformed = replace(current, coefficients=np.array([
            gauge.matrix(q).conj().T @ c for q, c in zip(current.q, current.coefficients)
        ]))
        identity = dict(family_id="analytic-singlet-triplet", trajectory_id="v304-soc", step=index)
        fixed = evaluate_molecular_spawn_checkpoint_v304(current, bridge, grid, **identity)
        moving = evaluate_molecular_spawn_checkpoint_v304(transformed, moving_bridge, moving_grid, **identity)
        refs, kinetic = _references(current, bridge, grid, fixed)
        mrefs, _ = _references(transformed, moving_bridge, moving_grid, moving)
        comparisons = []
        for a, b, ref, mref in zip(fixed["candidates"], moving["candidates"], refs, mrefs):
            ar, br = a["record"], b["record"]
            comparisons.append({
                "serial": a["serial"], "fixed_gain": ar["predicted_residual_reduction"],
                "moving_gain": br["predicted_residual_reduction"],
                "fixed_relative_gain": ar["relative_residual_reduction"],
                "moving_relative_gain": br["relative_residual_reduction"],
                "fixed_reference_gain": ref["complete_lstsq_gain"],
                "moving_reference_gain": mref["complete_lstsq_gain"],
                "fixed_legacy_gain": ref["legacy_real_only_gain"],
                "moving_legacy_gain": mref["legacy_real_only_gain"],
                "fixed_accepted": ar["accepted_by_residual_oracle"], "moving_accepted": br["accepted_by_residual_oracle"],
                "fixed_reasons": ar["rejection_reasons"], "moving_reasons": br["rejection_reasons"],
                "fixed_rank": ar["schur_rank"], "moving_rank": br["schur_rank"],
            })
        points.append({"time_au": float(current.time_au), "fixed": fixed, "moving": moving,
                       "kinetic_action_norm": kinetic, "comparisons": comparisons})
    replay = replay_molecular_spawn_checkpoint_v304(
        serialize_molecular_spawn_checkpoint_v304(points[-1]["moving"]), moving_bridge, moving_grid,
    )
    rejected = evaluate_molecular_spawn_checkpoint_v304(
        state, bridge, grid, family_id="analytic-singlet-triplet", trajectory_id="v304-soc", step=0,
        settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e6),
    )
    evidence = {
        "schema": SCHEMA_V304_EVIDENCE,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "snapshots": points, "propagation_receipt": step.receipt,
        "same_runtime_exact_replay": replay == points[-1]["moving"],
        "rejected_checkpoint": rejected,
        "dataset_certified": False, "live_molecular_soc_validated": False,
        "scope": "Local residual screening on one analytic SOC trajectory; fixed finite grid; no topology change, final admission, or certified training dataset.",
    }
    evidence["gates"] = evaluate_spawn_gates_v304(evidence)
    evidence["passed"] = all(evidence["gates"].values())
    evidence["fingerprint"] = stable_hash_v294(evidence)
    return evidence
