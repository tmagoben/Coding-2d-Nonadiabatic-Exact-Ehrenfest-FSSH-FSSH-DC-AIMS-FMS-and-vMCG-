"""Fresh unseeded birth, activation and repeated-admission lifecycle evidence."""
from dataclasses import replace
import json
import platform
import re

import numpy as np
import scipy
from scipy.sparse.linalg import expm_multiply

from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import evaluate_curved_state_v293
from .molecular_active_v306 import molecular_active_midpoint_step_v306
from .molecular_adaptive_v306 import (
    run_molecular_adaptive_v306, serialize_adaptive_trajectory_v306,
    replay_adaptive_trajectory_v306,
)
from .molecular_admission_v305 import (
    evaluate_molecular_admission_v305, serialize_molecular_admission_v305,
    replay_molecular_admission_v305,
)
from .molecular_spawn_v304 import _checkpoint_connection
from .pre_ml_contract_v294 import stable_hash_v294
from .v302_benchmark import free_packet_fixture_v302

SCHEMA_V306_EVIDENCE = "gnd-molecular-adaptive-evidence-v0.30.6-dev"
TIME_STEPS_V306 = (.02,)*12 + (.001,)*4 + (.004,)*4 + (.0001,)*2
SHAPES = ("q", "p", "width_matrices", "chirp_matrices")


def _sealed(receipt):
    payload = dict(receipt)
    digest = payload.pop("fingerprint", None)
    return stable_hash_v294(payload) == digest


def _coefficients(state):
    return np.array(state["coefficients_real"]) + 1j*np.array(state["coefficients_imag"])


def _field_payload(psi):
    return {"real": psi.real.tolist(), "imag": psi.imag.tolist()}


def adaptive_measurements_v306(evidence):
    """Derive reported errors and drifts from saved fields and step observables."""
    reference = evidence["reference"]
    def field(name):
        data = reference[name]
        return np.array(data["real"]) + 1j*np.array(data["imag"])
    def norm(x):
        return float(np.sqrt(np.sum(np.array(reference["weights"])[:, None] * abs(x)**2)))
    steps = evidence["trajectory"]["steps"]
    first = steps[0]["propagation"]["before_observables"]
    observables = [s["propagation"]["after_observables"] for s in steps]
    return {
        "adaptive_grid_exact_error": norm(field("adaptive") - field("exact")),
        "fixed_packet_grid_exact_error": norm(field("fixed_packet") - field("exact")),
        "maximum_norm_drift": max(abs(o["norm"]-first["norm"]) for o in observables),
        "maximum_energy_drift_hartree": max(abs(o["energy_hartree"]-first["energy_hartree"]) for o in observables),
        "maximum_scaled_midpoint_residual": max(s["propagation"]["scaled_midpoint_residual"] for s in steps),
        "final_time_au": steps[-1]["after"]["time_au"],
    }


def evaluate_adaptive_gates_v306(evidence):
    """Reconstruct acceptance from measurements, never from reported gate flags.

    This offline audit cannot replace the physical replay executed by the campaign.
    Hashes detect accidental edits; the replay reruns every physical decision.
    """
    trajectory = evidence["trajectory"]
    steps = trajectory["steps"]
    initial = evidence["initial_admission"]
    propagations = [s["propagation"] for s in steps]
    events = [initial] + [s["topology"] for s in steps if s["topology"] is not None]
    activations = [(s, row) for s in steps for row in s["activation"]["trials"] if row["activated"]]
    sealed = [trajectory, initial] + steps + propagations + [s["activation"] for s in steps] + events
    measured = adaptive_measurements_v306(evidence)

    def lifecycle_bound():
        prior_state, prior_metadata = initial["after_state"], initial["after_metadata"]
        for index, step in enumerate(steps, 1):
            propagation, activation = step["propagation"], step["activation"]
            if step["before"] != prior_state or step["before_metadata"] != prior_metadata:
                return False
            if propagation["before"] != step["before"] or propagation["dt"] != step["dt"]:
                return False
            if propagation["active_shape_mask"] != prior_metadata["active_shape_mask"] or activation["mask_before"] != prior_metadata["active_shape_mask"]:
                return False
            if activation["state_fingerprint"] != stable_hash_v294(propagation["after"]):
                return False
            expected = {**prior_metadata, "packet_ages": [age+1 for age in prior_metadata["packet_ages"]],
                        "active_shape_mask": activation["mask_after"]}
            topology = step["topology"]
            if topology is None:
                if step["after"] != propagation["after"] or step["after_metadata"] != expected:
                    return False
            elif (topology["screening"]["state"] != propagation["after"] or topology["before_metadata"] != expected
                  or step["after"] != topology["after_state"] or step["after_metadata"] != topology["after_metadata"]):
                return False
            if step["identity"] != {**initial["screening"]["identity"], "step": index}:
                return False
            if abs(propagation["after"]["time_au"] - step["before"]["time_au"] - step["dt"]) > 1e-14:
                return False
            prior_state, prior_metadata = step["after"], step["after_metadata"]
        return True

    def activation_passes(step, row):
        config = step["basis_settings"]
        spectrum = np.asarray(row["projected_shape_metric_spectrum"])
        if spectrum.ndim != 1 or spectrum.size == 0 or not np.all(np.isfinite(spectrum)):
            return False
        cutoff = max(config["tangent_absolute_cutoff"], config["tangent_relative_cutoff"]*spectrum[0])
        rank = int(np.count_nonzero(spectrum > cutoff))
        d = len(step["after"]["q"][0])
        speed = row["trial_shape_speed"] / max(row["base_shape_speed"], 1e-30)
        coefficient_norm = float(np.linalg.norm(_coefficients(step["propagation"]["after"])[row["packet_index"]])**2)
        return bool(rank == row["projected_shape_rank"] == 2*d+d*(d+1)
                    and np.isclose(coefficient_norm, row["coefficient_norm_squared"], rtol=1e-12, atol=1e-20)
                    and coefficient_norm >= config["inherited"]["newborn_activation_population"]
                    and row["projected_shape_cutoff"] == cutoff
                    and np.isclose(spectrum[0]/spectrum[-1], row["projected_shape_condition"])
                    and spectrum[0]/spectrum[-1] <= config["inherited"]["maximum_activation_condition_number"]
                    and np.isclose(speed, row["shape_velocity_amplification"])
                    and np.isfinite(speed) and speed <= config["inherited"]["maximum_activation_velocity_amplification"])

    def event_passes(event):
        selected = [r["attempt"] for r in event["candidate_decisions"] if r["final_admitted"]]
        config = event["screening"]["settings"]["inherited"]
        return (event["event_kind"] == "spawn" and len(selected) == 1 and all(
            row is not None and row["projection"] is not None
            and 0 <= row["projection"]["relative_loss"] <= config["maximum_projection_loss"]
            and abs(row["energy_after_hartree"]-row["energy_before_hartree"]) <= config["maximum_event_energy_jump_hartree"]
            and np.isclose(row["energy_jump_hartree"], row["energy_after_hartree"]-row["energy_before_hartree"], atol=1e-15, rtol=1e-10)
            for row in selected))

    moved = any(any(np.linalg.norm(np.array(later["propagation"]["after"][name])[row["packet_index"]]
                                  - np.array(step["after"][name])[row["packet_index"]]) > 1e-6 for name in SHAPES)
                for step, row in activations for later in steps if later["identity"]["step"] > step["identity"]["step"])
    gates = {
        "source_revision_recorded": re.fullmatch(r"[0-9a-f]{40}", evidence["generation_commit"]) is not None,
        "sealed_receipts": all(_sealed(s) for s in sealed),
        "complete_explicit_schedule": trajectory["time_steps"] == list(TIME_STEPS_V306) == [s["dt"] for s in steps],
        "state_and_metadata_continuity": lifecycle_bound(),
        "dormant_shapes_exactly_frozen": all(np.array_equal(np.array(p["before"][name])[~np.array(p["active_shape_mask"])], np.array(p["after"][name])[~np.array(p["active_shape_mask"])]) for p in propagations for name in SHAPES),
        "unseeded_newborns": all(np.linalg.norm(_coefficients(e["after_state"])[-1])**2 < 1e-20 and not e["after_metadata"]["active_shape_mask"][-1] and e["after_metadata"]["packet_ages"][-1] == 0 for e in events),
        "both_newborns_gain_amplitude": len(steps[-1]["after"]["q"]) == 3 and np.all(np.sum(abs(_coefficients(steps[-1]["after"])[1:])**2, axis=1) > 1e-10),
        "activation_uses_measured_gates": len(activations) >= 1 and all(activation_passes(s, r) for s, r in activations),
        "activated_shape_subsequently_moves": moved,
        "two_admitted_spawns": len(events) == 2 and all(event_passes(e) for e in events),
        "topology_deferred_for_dormancy": all((s["topology_status"] == "deferred_dormant_shapes" and s["topology"] is None) if not all(s["activation"]["mask_after"]) else s["topology_status"] == "evaluated" for s in steps),
        "accepted_nonlinear_and_metric_residuals": all(0 <= p["scaled_midpoint_residual"] <= 1 and 0 <= p["midpoint_svd_residual"] < 2e-8 and 0 <= p["midpoint_null_rhs_relative"] < 2e-8 for p in propagations),
        "bounded_norm_and_energy_drift": measured["maximum_norm_drift"] < 2e-7 and measured["maximum_energy_drift_hartree"] < 1e-7,
        "bounded_fixed_grid_reference_error": measured["adaptive_grid_exact_error"] < .002,
        "same_runtime_physical_replay": evidence["initial_admission_replayed"] is True and evidence["trajectory_replayed"] is True,
        "no_training_certification": not evidence["dataset_certified"] and not evidence["training_eligible"] and all(not s["dataset_certified"] and not s["training_eligible"] for s in steps),
    }
    return {key: bool(value) for key, value in gates.items()}


def run_v306_adaptive_campaign(*, generation_commit, verbose=False):
    if not isinstance(generation_commit, str) or re.fullmatch(r"[0-9a-f]{40}", generation_commit) is None:
        raise ValueError("generation_commit must be a full source commit SHA")
    state, bridge, grid = free_packet_fixture_v302()
    default = CurvedBasisSettingsV292()
    settings = replace(default, inherited=replace(default.inherited, maximum_packet_count=3, adapt_every_steps=1))
    identity = dict(family_id="analytic-free-grid", trajectory_id="v306-unseeded-lifecycle")
    initial = evaluate_molecular_admission_v305(state, bridge, grid, **identity, step=0, settings=settings)
    callback = (lambda row: print(json.dumps(row, sort_keys=True), flush=True)) if verbose else None
    trajectory = run_molecular_adaptive_v306(initial.after, bridge, grid, time_steps=TIME_STEPS_V306,
                                            **identity, **initial.receipt["after_metadata"],
                                            basis_settings=settings, progress=callback)
    serialized = serialize_adaptive_trajectory_v306(trajectory)
    if verbose:
        print("Replaying initial admission and every adaptive solve and decision", flush=True)
    replayed_initial = replay_molecular_admission_v305(serialize_molecular_admission_v305(initial.receipt), bridge, grid)
    replayed = replay_adaptive_trajectory_v306(serialized, bridge, grid)
    if verbose:
        print("Propagating a fixed one-packet baseline on the same time schedule", flush=True)
    baseline = state
    for dt in TIME_STEPS_V306:
        baseline = molecular_active_midpoint_step_v306(baseline, bridge, grid, dt).after
    connection = _checkpoint_connection(bridge)
    field = lambda s: evaluate_curved_state_v293(s, connection, grid.quadrature)
    exact = expm_multiply(-1j*trajectory[-1].after.time_au*grid.matrix, field(state).ravel()).reshape(-1, state.nstate)
    evidence = {
        "schema": SCHEMA_V306_EVIDENCE, "generation_commit": generation_commit,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "operator": grid.evidence(), "initial_admission": initial.receipt, "trajectory": json.loads(serialized),
        "initial_admission_replayed": replayed_initial.receipt == initial.receipt,
        "trajectory_replayed": serialize_adaptive_trajectory_v306(replayed) == serialized,
        "reference": {"weights": grid.quadrature.weights.tolist(), "exact": _field_payload(exact),
                      "adaptive": _field_payload(field(trajectory[-1].after)), "fixed_packet": _field_payload(field(baseline)),
                      "time_steps": list(TIME_STEPS_V306)},
        "dataset_certified": False, "training_eligible": False,
        "scope": "One scalar analytic finite-grid lifecycle; explicit smaller post-activation steps; no adaptive time-step controller, convergence certification, live molecular SOC, or training dataset.",
    }
    evidence["measurements"] = adaptive_measurements_v306(evidence)
    evidence["gates"] = evaluate_adaptive_gates_v306(evidence)
    evidence["passed"] = all(evidence["gates"].values())
    evidence["fingerprint"] = stable_hash_v294(evidence)
    return evidence
