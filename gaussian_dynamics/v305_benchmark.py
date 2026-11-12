"""Fresh analytic SOC admission and replay-verified export campaign."""
from dataclasses import replace
import platform

import numpy as np
import scipy

from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import evaluate_curved_state_v293
from .molecular_admission_v305 import evaluate_molecular_admission_v305, serialize_molecular_admission_v305
from .molecular_dataset_v305 import MolecularAdmissionSourceV305, materialize_admission_dataset_v305
from .molecular_spawn_v304 import _checkpoint_connection
from .pre_ml_contract_v294 import stable_hash_v294
from .v303_benchmark import soc_packet_fixture_v303

SCHEMA_V305_EVIDENCE = "gnd-molecular-admission-evidence-v0.30.5-dev"


def _sealed(receipt):
    payload = dict(receipt)
    digest = payload.pop("fingerprint", None)
    return stable_hash_v294(payload) == digest


def evaluate_admission_gates_v305(evidence):
    points = evidence["cases"]
    receipts = [point[frame] for point in points for frame in ("fixed", "moving")]
    dataset = evidence["dataset"]
    rows = dataset["rows"]
    def near(a, b):
        return np.isfinite(a) and np.isfinite(b) and abs(a-b) <= 3e-10 + 3e-6 * max(abs(a), abs(b))
    def decisions_agree(point):
        return all(a["final_admitted"] == b["final_admitted"] and a["outcome"] == b["outcome"]
                   for a, b in zip(point["fixed"]["candidate_decisions"], point["moving"]["candidate_decisions"]))
    def projections_pass(receipt):
        attempts = [row["attempt"] for row in receipt["candidate_decisions"] if row["final_admitted"]]
        controls = receipt["screening"]["settings"]["inherited"]
        return len(attempts) == 1 and all(
            a is not None and a["projection"] is not None
            and 0 <= a["projection"]["relative_loss"] <= controls["maximum_projection_loss"]
            and abs(a["energy_jump_hartree"]) <= controls["maximum_event_energy_jump_hartree"] for a in attempts)
    def features_agree(point):
        identity = point["fixed"]["screening"]["identity"]
        variants = {frame: [row for row in rows if row["identity"]["trajectory_id"] == identity["trajectory_id"] and row["identity"]["variant_id"] == frame] for frame in ("fixed", "moving")}
        if len(variants["fixed"]) != 4 or len(variants["moving"]) != 4:
            return False
        for fixed, moving in zip(variants["fixed"], variants["moving"]):
            if fixed["features"].keys() != moving["features"].keys():
                return False
            for key, a in fixed["features"].items():
                b = moving["features"][key]
                if isinstance(a, str):
                    if a != b:
                        return False
                elif not np.allclose(a, b, atol=3e-10, rtol=3e-6):
                    return False
        return True
    controls = evidence["controls"]
    by_fingerprint = {r["fingerprint"]: r for r in receipts + list(controls.values())}
    def labels_bound():
        if sum(r["screening"]["candidate_count"] for r in by_fingerprint.values()) != len(rows):
            return False
        for row in rows:
            receipt = by_fingerprint.get(row["provenance"]["admission_fingerprint"])
            if receipt is None:
                return False
            serial = row["identity"]["candidate_serial"]
            if serial >= len(receipt["candidate_decisions"]):
                return False
            decision = receipt["candidate_decisions"][serial]
            native = receipt["screening"]["candidates"][serial]["record"]
            if any(row["targets"][key] != decision[key] for key in ("final_admitted", "outcome", "residual_accepted")):
                return False
            if row["targets"]["predicted_residual_reduction"] != native["predicted_residual_reduction"]:
                return False
        return True
    groups = {}
    for row in rows:
        key = (row["identity"]["family_id"], row["identity"]["trajectory_id"])
        groups.setdefault(key, set()).add(row["split"])
    gates = {
        "two_distinct_initial_conditions": len(points) == 2 and points[0]["fixed"]["screening"]["state"]["q"] != points[1]["fixed"]["screening"]["state"]["q"],
        "sealed_evidence": all(_sealed(r) for r in by_fingerprint.values()) and _sealed(dataset) and all(_sealed(row) for row in rows),
        "all_candidates_recorded": all(len(r["candidate_decisions"]) == r["screening"]["candidate_count"] == 4 for r in receipts),
        "exactly_one_spawn": all(r["event_kind"] == "spawn" and sum(row["final_admitted"] for row in r["candidate_decisions"]) == 1 for r in receipts),
        "projection_and_energy_gates": all(projections_pass(r) for r in receipts),
        "newborn_shapes_dormant": all(r["after_metadata"]["active_shape_mask"] == [True, False] and r["after_metadata"]["packet_ages"] == [0, 0] for r in receipts),
        "physical_gauge_scores": all(near(a["record"]["predicted_residual_reduction"], b["record"]["predicted_residual_reduction"]) for p in points for a, b in zip(p["fixed"]["screening"]["candidates"], p["moving"]["screening"]["candidates"])),
        "physical_gauge_admissions": all(decisions_agree(p) for p in points),
        "physical_gauge_after_state": all(np.isfinite(p["after_wavefunction_gauge_error"]) and 0 <= p["after_wavefunction_gauge_error"] < 3e-9 for p in points),
        "physical_gauge_export_features": all(features_agree(p) for p in points),
        "projection_rejection_retained": any(row["outcome"] == "projection_rejected" and not row["final_admitted"] for row in controls["projection"]["candidate_decisions"]),
        "lifecycle_block_retained": all(row["outcome"] == "maximum_packet_count" and not row["final_admitted"] for row in controls["packet_limit"]["candidate_decisions"]),
        "residual_rejection_retained": all(row["outcome"] == "residual_rejected" and not row["final_admitted"] for row in controls["residual"]["candidate_decisions"]),
        "all_checkpoints_physically_replayed": dataset["exact_replay_verified"] is True and dataset["checkpoint_count"] == len(by_fingerprint) == 7,
        "export_labels_bound_to_admission": labels_bound(),
        "export_manifest_counts": dataset["record_count"] == len(rows) == 28 and dataset["checkpoint_count"] == len(dataset["checkpoints"]) == 7 and sum(dataset["split_counts"].values()) == len(rows),
        "trajectory_splits_do_not_leak": len(groups) == 2 and all(len(splits) == 1 for splits in groups.values()),
        "no_invented_target_or_certification": all(r["target_state"] is None and not r["dataset_certified"] for r in receipts) and all(row["target_state"] is None for row in rows) and not any(dataset[k] for k in ("dataset_certified", "physical_invariance_certified", "training_eligible")),
    }
    return {key: bool(value) for key, value in gates.items()}


def run_v305_admission_campaign(*, generation_commit, verbose=False):
    original, bridge, grid, gauge, _ = soc_packet_fixture_v303()
    _, moving_bridge, moving_grid, _, _ = soc_packet_fixture_v303(moving=True)
    points, sources = [], []
    for index, state in enumerate((original, replace(original, q=np.array([[-.25]]), p=np.array([[.55]])))):
        if verbose:
            print(f"Admission and gauge audit for initial condition {index}", flush=True)
        identity = dict(family_id="analytic-singlet-triplet", trajectory_id=f"v305-initial-{index}", step=0)
        transformed = replace(state, coefficients=np.array([gauge.matrix(q).conj().T @ c for q, c in zip(state.q, state.coefficients)]))
        fixed = evaluate_molecular_admission_v305(state, bridge, grid, **identity)
        moving = evaluate_molecular_admission_v305(transformed, moving_bridge, moving_grid, **identity)
        a = evaluate_curved_state_v293(fixed.after, _checkpoint_connection(bridge), grid.quadrature)
        b = evaluate_curved_state_v293(moving.after, _checkpoint_connection(moving_bridge), moving_grid.quadrature)
        b = np.array([gauge.matrix(q) @ psi for q, psi in zip(grid.quadrature.points, b)])
        error = float(np.sqrt(np.sum(grid.quadrature.weights[:, None] * abs(a-b)**2)))
        points.append({"fixed": fixed.receipt, "moving": moving.receipt, "after_wavefunction_gauge_error": error})
        sources.extend((MolecularAdmissionSourceV305(serialize_molecular_admission_v305(fixed.receipt), bridge, grid, "fixed"),
                        MolecularAdmissionSourceV305(serialize_molecular_admission_v305(moving.receipt), moving_bridge, moving_grid, "moving")))
    default = CurvedBasisSettingsV292()
    settings = {"projection": replace(default, inherited=replace(default.inherited, maximum_basis_condition_number=1.01)),
                "packet_limit": replace(default, inherited=replace(default.inherited, maximum_packet_count=1)),
                "residual": replace(default, minimum_predicted_residual_reduction=1e6)}
    controls = {}
    for kind, config in settings.items():
        if verbose:
            print(f"Evaluating {kind} rejection control", flush=True)
        result = evaluate_molecular_admission_v305(original, bridge, grid, **points[0]["fixed"]["screening"]["identity"], settings=config)
        controls[kind] = result.receipt
        sources.append(MolecularAdmissionSourceV305(serialize_molecular_admission_v305(result.receipt), bridge, grid, kind))
    if verbose:
        print("Replaying every checkpoint before dataset export", flush=True)
    dataset = materialize_admission_dataset_v305(sources, generation_commit=generation_commit)
    evidence = {"schema": SCHEMA_V305_EVIDENCE, "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
                "cases": points, "controls": controls, "dataset": dataset,
                "scope": "Two analytic SOC initial conditions on a fixed finite grid; one topology event per checkpoint; no adaptive propagation or training certification."}
    evidence["gates"] = evaluate_admission_gates_v305(evidence)
    evidence["passed"] = all(evidence["gates"].values())
    evidence["fingerprint"] = stable_hash_v294(evidence)
    return evidence
