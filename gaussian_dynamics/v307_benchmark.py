"""Adaptive SOC time refinement, gauge convergence and checkpoint sensitivity.

The stored evidence contains checkpoint states, every microstep's diagnostics,
physical fields and full-trace digests. Selected complete traces are physically
replayed during generation; offline gates do not pretend to rerun that replay.
"""
from dataclasses import replace
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from pathlib import Path
import json
import platform
import re

import numpy as np
import scipy
from scipy.sparse.linalg import expm_multiply

from .analytic_soc_bridge_v307 import analytic_soc_fixture_v307
from .curved_basis_adaptation_v292 import CurvedBasisSettingsV292
from .curved_basis_adaptation_v293 import evaluate_curved_state_v293
from .molecular_admission_v305 import evaluate_molecular_admission_v305, serialize_molecular_admission_v305, replay_molecular_admission_v305
from .molecular_checkpoints_v307 import run_molecular_checkpoints_v307, serialize_checkpoint_trajectory_v307, replay_checkpoint_trajectory_v307
from .pre_ml_contract_v294 import stable_hash_v294
from .v306_benchmark import _field_payload, _coefficients, _sealed, SHAPES

SCHEMA_V307_EVIDENCE = "gnd-adaptive-soc-convergence-v0.30.7-dev"
CHECKPOINTS_V307 = (.02, .04, .06, .08)
SUBSTEPS_V307 = (20, 40, 80)
CASES_V307 = (("central", .1, .3), ("displaced", -.25, .55))
CADENCES_V307 = (.02, .01, .005, .0025)


def _field(payload):
    return np.array(payload["real"]) + 1j*np.array(payload["imag"])


def _norm(field, weights):
    return float(np.sqrt(np.sum(np.asarray(weights)[:, None]*abs(field)**2)))


def _triplet(field, weights):
    return float(np.sum(np.asarray(weights)[:, None]*abs(field[:, 1:])**2) / _norm(field, weights)**2)


def _props(run):
    return [p for checkpoint in run["checkpoints"] for p in checkpoint["propagations"]]


def _activation_times(run):
    return [r["checkpoint_time"] for r in run["checkpoints"] if any(x["activated"] for x in r["activation"]["trials"])]


def _bracket(run):
    previous = 0.
    for row in run["checkpoints"]:
        if any(x["activated"] for x in row["activation"]["trials"]):
            return [previous, row["checkpoint_time"]]
        previous = row["checkpoint_time"]
    return None


def _run_metrics(run):
    props = _props(run)
    initial = props[0]["before_observables"]
    weights = run["weights"]
    final, exact = _field(run["final_field"]), _field(run["exact_field"])
    return {
        "grid_exact_error": _norm(final-exact, weights),
        "triplet_population": _triplet(final, weights),
        "exact_triplet_population": _triplet(exact, weights),
        "maximum_norm_drift": max(abs(p["after_observables"]["norm"]-initial["norm"]) for p in props),
        "maximum_energy_drift_hartree": max(abs(p["after_observables"]["energy_hartree"]-initial["energy_hartree"]) for p in props),
        "activation_times": _activation_times(run),
        "activation_bracket": _bracket(run),
        "maximum_scaled_midpoint_residual": max(p["scaled_midpoint_residual"] for p in props),
    }


def _order(differences):
    return float(np.log2(differences[0]/differences[1])) if all(d > 0 for d in differences) else None


def convergence_measurements_v307(evidence):
    cases = []
    for case in evidence["cases"]:
        frames = {}
        for name, frame in case["frames"].items():
            runs = frame["levels"]
            differences = [_norm(_field(a["final_field"])-_field(b["final_field"]), a["weights"]) for a,b in zip(runs[:-1], runs[1:])]
            frames[name] = {"levels": [_run_metrics(r) for r in runs], "endpoint_differences": differences, "observed_order": _order(differences)}
        fixed, moving = (case["frames"][n]["levels"] for n in ("fixed", "moving"))
        gauge_errors = [_norm(_field(a["final_field"])-_field(b["final_field"]), a["weights"]) for a,b in zip(fixed,moving)]
        cases.append({"name": case["name"], "frames": frames, "gauge_errors": gauge_errors,
                      "gauge_orders": [_order(gauge_errors[i:i+2]) for i in (0,1)]})
    cadence = evidence["cadence"]
    differences = [_norm(_field(a["final_field"])-_field(b["final_field"]), a["weights"]) for a,b in zip(cadence[:-1],cadence[1:])]
    all_runs = [r for c in evidence["cases"] for f in c["frames"].values() for r in f["levels"]] + cadence + [evidence["soc_off"]]
    return {"cases": cases, "cadence": {"spacings": list(CADENCES_V307), "levels": [_run_metrics(r) for r in cadence],
                                        "successive_endpoint_differences": differences},
            "soc_off": _run_metrics(evidence["soc_off"]),
            "maximum_scaled_midpoint_residual": max(_run_metrics(r)["maximum_scaled_midpoint_residual"] for r in all_runs)}


def _activation_valid(row, settings, state):
    cfg = settings["inherited"]
    spectrum = np.asarray(row["projected_shape_metric_spectrum"])
    cutoff = max(settings["tangent_absolute_cutoff"], settings["tangent_relative_cutoff"]*float(spectrum[0]))
    population = float(np.linalg.norm(_coefficients(state)[row["packet_index"]])**2)
    speed = row["trial_shape_speed"]/max(row["base_shape_speed"],1e-30)
    return (len(spectrum) == 4 and np.all(spectrum > cutoff) and row["projected_shape_rank"] == 4
            and np.isclose(row["coefficient_norm_squared"],population,atol=1e-20,rtol=1e-10)
            and population >= cfg["newborn_activation_population"]
            and row["projected_shape_cutoff"] == cutoff
            and np.isclose(row["projected_shape_condition"],spectrum[0]/spectrum[-1])
            and spectrum[0]/spectrum[-1] <= cfg["maximum_activation_condition_number"]
            and np.isclose(row["shape_velocity_amplification"],speed)
            and speed <= cfg["maximum_activation_velocity_amplification"])


def _schedule_bound(run, initial):
    prior, meta = initial["after_state"], initial["after_metadata"]
    times = run["checkpoint_times"]
    if len(run["checkpoints"]) != len(times):
        return False
    for index, (end, row) in enumerate(zip(times,run["checkpoints"]),1):
        if row["checkpoint_time"] != end or row["identity"]["step"] != index or row["before"] != prior or row["before_metadata"] != meta:
            return False
        if len(row["propagations"]) != run["substeps"] or row["substeps"] != run["substeps"]:
            return False
        targets = np.linspace(prior["time_au"],end,run["substeps"]+1)
        time = prior["time_au"]
        for target, p in zip(targets[1:],row["propagations"]):
            if (p["before_time"] != time or p["dt"] != float(target-time)
                or p["after_time"] != float(time+p["dt"]) or p["active_shape_mask"] != meta["active_shape_mask"]):
                return False
            time = p["after_time"]
        if time != end or row["after"]["time_au"] != end:
            return False
        activation = row["activation"]
        expected = {**meta, "packet_ages": [age+1 for age in meta["packet_ages"]], "active_shape_mask": activation["mask_after"]}
        if (activation["mask_before"] != meta["active_shape_mask"] or expected != row["after_metadata"]
            or activation["state_fingerprint"] != stable_hash_v294(row["after"])):
            return False
        if all(expected["active_shape_mask"]):
            if row["topology_status"] != "evaluated" or row["topology_event"] != "none" or row["topology_reason"] != "maximum_packet_count":
                return False
        elif row["topology_status"] != "deferred_dormant_shapes" or row["topology_event"] is not None:
            return False
        prior, meta = row["after"], row["after_metadata"]
    return True


def evaluate_convergence_gates_v307(evidence):
    """Derive bounded acceptance from measurements, never saved success flags."""
    cases, cadence, off = evidence["cases"], evidence["cadence"], evidence["soc_off"]
    frames = [f for c in cases for f in c["frames"].values()]
    pairs = [(r,f["initial_admission"]) for f in frames for r in f["levels"]]
    pairs += [(r,cases[0]["frames"]["fixed"]["initial_admission"]) for r in cadence]
    pairs += [(off,evidence["soc_off_initial_admission"])]
    runs = [r for r,_ in pairs]
    checkpoints = [c for r in runs for c in r["checkpoints"]]
    propagations = [p for r in runs for p in _props(r)]
    initial = [f["initial_admission"] for f in frames] + [evidence["soc_off_initial_admission"]]
    activations = [(c,x) for c in checkpoints for x in c["activation"]["trials"] if x["activated"]]
    measured = convergence_measurements_v307(evidence)
    def admission_passes(event):
        selected = [r["attempt"] for r in event["candidate_decisions"] if r["final_admitted"]]
        cfg = event["screening"]["settings"]["inherited"]
        return (event["event_kind"] == "spawn" and len(selected) == 1
                and np.linalg.norm(_coefficients(event["after_state"])[-1])**2 < 1e-20
                and event["after_metadata"]["active_shape_mask"] == [True,False]
                and all(x["projection"]["relative_loss"] <= cfg["maximum_projection_loss"]
                        and abs(x["energy_after_hartree"]-x["energy_before_hartree"]) <= cfg["maximum_event_energy_jump_hartree"] for x in selected))
    def moved(run):
        first = next((c for c in run["checkpoints"] if any(x["activated"] for x in c["activation"]["trials"])),None)
        return first is not None and any(np.linalg.norm(np.array(run["checkpoints"][-1]["after"][key])[1]-np.array(first["after"][key])[1]) > 1e-6 for key in SHAPES)
    fixed_runs = [r for f in frames for r in f["levels"]]
    brackets = [_bracket(r) for r in cadence]
    gates = {
        "source_revision_recorded": isinstance(evidence["generation_commit"],str) and re.fullmatch(r"[0-9a-f]{40}",evidence["generation_commit"]) is not None,
        "complete_two_case_two_gauge_campaign": [(c["name"],c["q0"],c["p0"]) for c in cases] == list(CASES_V307)
            and all(set(c["frames"]) == {"fixed","moving"} for c in cases)
            and all(tuple(r["substeps"] for r in f["levels"]) == SUBSTEPS_V307 for f in frames),
        "fixed_physical_checkpoint_schedule": all(tuple(r["checkpoint_times"]) == CHECKPOINTS_V307 for r in fixed_runs),
        "state_clock_and_metadata_continuity": all(_schedule_bound(r,i) for r,i in pairs),
        "sealed_admissions_and_activations": all(_sealed(r) for r in initial) and all(_sealed(c["activation"]) for c in checkpoints),
        "unseeded_admitted_births": all(admission_passes(i) for i in initial),
        "dormant_shapes_exactly_frozen": all(np.array_equal(np.array(c["before"][k])[~np.array(c["before_metadata"]["active_shape_mask"])],np.array(c["after"][k])[~np.array(c["before_metadata"]["active_shape_mask"])]) for c in checkpoints for k in SHAPES),
        "measured_activation_gates": all(len(_activation_times(r)) == 1 for r in fixed_runs+cadence) and all(_activation_valid(x,c["activation"]["settings"],c["after"]) for c,x in activations),
        "activated_shapes_continue_to_move": all(moved(r) for r in fixed_runs+cadence),
        "accepted_midpoint_metric_and_null_residuals": all(0 <= p["scaled_midpoint_residual"] <= 1 and 0 <= p["midpoint_svd_residual"] <= 2e-8 and 0 <= p["midpoint_null_rhs_relative"] <= 2e-8 and p["midpoint_metric_condition"] <= 1e10 for p in propagations),
        "bounded_norm_and_energy_drift": all(_run_metrics(r)["maximum_norm_drift"] < 2e-6 and _run_metrics(r)["maximum_energy_drift_hartree"] < 1e-6 for r in runs),
        "second_order_fixed_checkpoint_convergence": all(f["observed_order"] is not None and 1.8 < f["observed_order"] < 2.2 and f["endpoint_differences"][-1] > 1e-9 for c in measured["cases"] for f in c["frames"].values()),
        "activation_time_stable_under_mesh_refinement": all(_activation_times(r) == [.04] for r in fixed_runs),
        "physical_gauge_error_converges": all(c["gauge_errors"][0] < 2e-6 and c["gauge_errors"][-1] < 1e-7 and all(o is not None and 1.7 < o < 2.3 for o in c["gauge_orders"]) for c in measured["cases"]),
        "bounded_finite_grid_reference_error": all(_run_metrics(r)["grid_exact_error"] < .002 for r in runs),
        "physical_soc_transfer_and_disabled_control": all(1e-5 < _run_metrics(r)["triplet_population"] < .01 for r in fixed_runs)
            and _run_metrics(off)["triplet_population"] < 1e-20 and _run_metrics(off)["exact_triplet_population"] < 1e-20
            and not _activation_times(off) and all(c["after_metadata"]["active_shape_mask"] == [True,False] for c in off["checkpoints"]),
        "checkpoint_cadence_has_controlled_microsteps": len(cadence) == 4 and all(np.allclose(r["checkpoint_times"],np.linspace(spacing,.08,round(.08/spacing)),rtol=0,atol=1e-15) and r["substeps"] == round(spacing/.0005) and all(abs(p["dt"]-.0005) < 1e-15 for p in _props(r)) for spacing,r in zip(CADENCES_V307,cadence)),
        "activation_bracket_narrows_with_checkpoint_cadence": all(b is not None for b in brackets) and all(b[0]-1e-14 <= c[0] < c[1] <= b[1]+1e-14 for b,c in zip(brackets[:-1],brackets[1:])) and np.isclose(brackets[-1][1]-brackets[-1][0],.0025),
        "bounded_checkpoint_policy_sensitivity": max(measured["cadence"]["successive_endpoint_differences"]) < .001,
        "selected_complete_physical_replays": all(f["initial_admission_replayed"] and f["levels"][0]["full_replay_matches"] is True for f in frames) and off["full_replay_matches"] is True and evidence["soc_off_initial_replayed"] is True,
        "no_training_or_continuum_certification": evidence["dataset_certified"] is False and evidence["training_eligible"] is False and evidence["continuum_certified"] is False and all(c["dataset_certified"] is False and c["training_eligible"] is False for c in checkpoints),
    }
    return {key: bool(value) for key,value in gates.items()}


def _settings_v307():
    default = CurvedBasisSettingsV292()
    return replace(default,inherited=replace(default.inherited,maximum_packet_count=2,adapt_every_steps=1))


def _setup_v307(q,p,moving=False,soc=True):
    state,bridge,grid,gauge,_ = analytic_soc_fixture_v307(moving=moving,soc_enabled=soc,q0=q,p0=p)
    identity = dict(family_id="analytic-singlet-triplet",trajectory_id=f"v307-{q}-{p}-{'moving' if moving else 'fixed'}-soc{soc}")
    initial = evaluate_molecular_admission_v305(state,bridge,grid,**identity,step=0,settings=_settings_v307())
    replay = replay_molecular_admission_v305(serialize_molecular_admission_v305(initial.receipt),bridge,grid)
    return state,bridge,grid,gauge,identity,initial,replay.receipt == initial.receipt

def _evolve_v307(model,times,n,*,replay=False,verbose=False):
    state,bridge,grid,gauge,identity,initial,_ = model
    trace = run_molecular_checkpoints_v307(initial.after,bridge,grid,checkpoint_times=times,substeps=n,**identity,
                                          **initial.receipt["after_metadata"],basis_settings=_settings_v307())
    serialized = serialize_checkpoint_trajectory_v307(trace)
    replay_match = None
    if replay:
        replay_match = serialize_checkpoint_trajectory_v307(replay_checkpoint_trajectory_v307(serialized,bridge,grid)) == serialized
    connection = bridge.connection_field()
    def physical(values):
        return np.array([gauge.matrix(q)@v for q,v in zip(grid.quadrature.points,values)]) if bridge.gauge is not None else values
    psi = evaluate_curved_state_v293(state,connection,grid.quadrature)
    exact = expm_multiply(-1j*times[-1]*grid.matrix,psi.ravel()).reshape(-1,4)
    rows = []
    for result in trace:
        r = result.receipt
        row = {k:r[k] for k in ("identity","checkpoint_time","substeps","before","after","before_metadata","after_metadata","activation","topology_status","dataset_certified","training_eligible")}
        row.update(topology_event=r["topology"]["event_kind"] if r["topology"] else None,
                   topology_reason=r["topology"]["event_reason"] if r["topology"] else None)
        row["propagations"] = [{**{k:p[k] for k in ("dt","iterations","active_shape_mask","scaled_midpoint_residual","midpoint_metric_condition","midpoint_metric_rank","midpoint_svd_residual","midpoint_null_rhs_relative","before_observables","after_observables")},
                                "before_time":p["before"]["time_au"],"after_time":p["after"]["time_au"]} for p in r["propagations"]]
        rows.append(row)
    run = {"checkpoint_times":list(times),"substeps":n,"checkpoints":rows,"full_trace_fingerprint":json.loads(serialized)["fingerprint"],"full_replay_matches":replay_match,
           "weights":grid.quadrature.weights.tolist(),"initial_field":_field_payload(physical(psi)),
           "final_field":_field_payload(physical(evaluate_curved_state_v293(trace[-1].after,connection,grid.quadrature))),"exact_field":_field_payload(physical(exact))}
    if verbose:
        print(json.dumps({"trajectory":identity["trajectory_id"],"substeps":n,"checkpoints":len(times),**_run_metrics(run)},sort_keys=True),flush=True)
    return run

def _frame_v307(arguments):
    name,q,p,moving,verbose = arguments
    model = _setup_v307(q,p,moving)
    return {"operator":model[2].evidence(),"initial_admission":model[5].receipt,"initial_admission_replayed":model[6],
            "levels":[_evolve_v307(model,CHECKPOINTS_V307,n,replay=n==SUBSTEPS_V307[0],verbose=verbose) for n in SUBSTEPS_V307]}


def _control_v307(arguments):
    spacing,verbose = arguments
    if spacing is None:
        model = _setup_v307(.1,.3,soc=False)
        return {"run":_evolve_v307(model,CHECKPOINTS_V307,20,replay=True,verbose=verbose),
                "initial_admission":model[5].receipt,"initial_replayed":model[6]}
    model = _setup_v307(.1,.3)
    return _evolve_v307(model,np.linspace(spacing,.08,round(.08/spacing)).tolist(),round(spacing/.0005),verbose=verbose)


def _campaign_part_v307(task):
    kind, arguments, revision, directory = task
    key = {"schema":SCHEMA_V307_EVIDENCE,"kind":kind,"arguments":arguments,"generation_commit":revision,
           "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    key = json.loads(json.dumps(key))
    path = Path(directory)/(stable_hash_v294(key)+".json") if directory is not None else None
    if path is not None and path.exists():
        payload = json.loads(path.read_text())
        if not _sealed(payload) or payload.get("key") != key:
            raise ValueError("campaign checkpoint integrity or provenance mismatch")
        return payload["result"]
    result = (_frame_v307 if kind == "frame" else _control_v307)(arguments)
    if path is not None:
        payload = {"key":key,"result":result}
        payload["fingerprint"] = stable_hash_v294(payload)
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(payload,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n")
        temporary.replace(path)
    return result


def run_v307_convergence_campaign(*, generation_commit, verbose=False, workers=2, checkpoint_directory=None):
    if not isinstance(generation_commit,str) or re.fullmatch(r"[0-9a-f]{40}",generation_commit) is None:
        raise ValueError("generation_commit must be a full source commit SHA")
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    if checkpoint_directory is not None:
        checkpoint_directory = str(Path(checkpoint_directory).resolve())
        Path(checkpoint_directory).mkdir(parents=True,exist_ok=True)
    frame_args = [(name,q,p,moving,verbose) for name,q,p in CASES_V307 for moving in (False,True)]
    control_args = [(spacing,verbose) for spacing in CADENCES_V307[1:]] + [(None,verbose)]
    frame_tasks = [("frame",args,generation_commit,checkpoint_directory) for args in frame_args]
    control_tasks = [("control",args,generation_commit,checkpoint_directory) for args in control_args]
    if workers == 1:
        frames = list(map(_campaign_part_v307,frame_tasks))
        controls = list(map(_campaign_part_v307,control_tasks))
    else:
        # Independent complete trajectories only; no parallel reductions inside a solve.
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
            frames = list(pool.map(_campaign_part_v307,frame_tasks))
            controls = list(pool.map(_campaign_part_v307,control_tasks))
    cases = [{"name":name,"q0":q,"p0":p,"frames":{"fixed":frames[2*i],"moving":frames[2*i+1]}}
             for i,(name,q,p) in enumerate(CASES_V307)]
    cadence = [cases[0]["frames"]["fixed"]["levels"][1],*controls[:-1]]
    off = controls[-1]
    evidence = {"schema":SCHEMA_V307_EVIDENCE,"generation_commit":generation_commit,
                "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
                "cases":cases,"cadence":cadence,"soc_off":off["run"],
                "soc_off_initial_admission":off["initial_admission"],"soc_off_initial_replayed":off["initial_replayed"],
                "dataset_certified":False,"training_eligible":False,"continuum_certified":False,
                "scope":"Two initial conditions, four-state 1D analytic SOC, two-packet budget, 13-point finite grid, t=0.08. Fixed-checkpoint time refinement is separate from sampled activation-time sensitivity. No live molecular provider, continuum convergence, general nonlinear solver, or training certification."}
    evidence["measurements"] = convergence_measurements_v307(evidence)
    evidence["gates"] = evaluate_convergence_gates_v307(evidence)
    evidence["passed"] = all(evidence["gates"].values())
    evidence["fingerprint"] = stable_hash_v294(evidence)
    return evidence
