"""Analytic SOC trajectory validation on a fixed one-dimensional finite grid.

Temporal self-convergence, representation covariance and Gaussian projection
error are measured separately. No continuum or live-molecular claim is made.
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import json
import platform

import numpy as np
import scipy
from scipy.sparse.linalg import expm_multiply

from .analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220, SingletTripletSOCConfigV220
from .complex_gauge_v21 import GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21, random_unitary_v21
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .molecular_frame_bridge_v300 import MolecularFrameBridgeV300
from .molecular_hamiltonian_v301 import build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301
from .molecular_midpoint_v302 import (
    run_molecular_midpoint_v302, serialize_molecular_trajectory_v302,
    replay_molecular_trajectory_v302,
)

SCHEMA_V303 = "gnd-soc-trajectory-evidence-v0.30.3-dev"


def soc_packet_fixture_v303(*, moving=False, soc_enabled=True):
    """A pure singlet with deliberately accelerated, analytic model parameters."""
    config = SingletTripletSOCConfigV220(
        mass=2.0, soc_scale=0.2 if soc_enabled else 0.0, soc_enabled=soc_enabled,
    )
    provider = AnalyticSingletTripletSOCProviderV220(config)
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, seed=303),
        np.array([[0.3], [-0.2], [0.1], [0.25]]), np.zeros(4),
    )
    bridge = MolecularFrameBridgeV300(
        GaugeTransformedOperatorProviderV21(provider, gauge) if moving else provider,
        [0.0],
    )
    grid = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-3, 3, 13),))
    state = CorrelatedGaussianSpinorStateV270(
        [[0.1]], [[0.3]], [[[1.1]]], [[[0.08]]], [[1.0, 0.0, 0.0, 0.0]],
    ).normalized()
    if moving:
        state = replace(state, coefficients=np.array([
            gauge.matrix(state.q[0]).conj().T @ state.coefficients[0],
        ]))
    return state, bridge, grid, gauge, config


def _weighted_norm(psi, weights):
    return float(np.sqrt(np.sum(weights[:, None] * abs(psi)**2)))


def _norm(psi, grid):
    return _weighted_norm(psi, grid.quadrature.weights)


def _physical_wavefunction(psi, grid, gauge, moving):
    # Transform the state into the fixed spin-diabatic basis. Raw component
    # populations in a mixed electronic gauge are not physical spin populations.
    if not moving:
        return psi
    return np.array([gauge.matrix(q) @ value for q, value in zip(grid.quadrature.points, psi)])


def _triplet_population(psi, grid):
    return float(np.sum(grid.quadrature.weights[:, None] * abs(psi[:, 1:])**2) / _norm(psi, grid)**2)


def _order(differences):
    if not all(np.isfinite(x) and x > 0 for x in differences):
        return None
    return float(np.log2(differences[0] / differences[1]))


def evaluate_soc_gates_v303(evidence):
    """Recompute acceptance from measurements, never from stored gate flags."""
    fixed, moving = (evidence["frames"][key] for key in ("fixed", "moving"))
    rows = fixed["refinement"] + moving["refinement"]
    finite = lambda value: value is not None and np.isfinite(value)
    below = lambda value, limit: finite(value) and 0 <= value < limit
    order_ok = lambda value: finite(value) and 1.8 < value < 2.2
    def frame_order_ok(frame):
        measured = _order(frame["successive_wavefunction_differences"])
        reported = frame["observed_time_order"]
        return order_ok(measured) and order_ok(reported) and abs(measured-reported) < 1e-12
    gauges = evidence["gauge_wavefunction_errors"]
    off = evidence["soc_disabled"]
    return {key: bool(value) for key, value in {
        "fixed_frame_second_order": frame_order_ok(fixed),
        "moving_frame_second_order": frame_order_ok(moving),
        "accepted_midpoint_residuals": all(finite(row["maximum_scaled_midpoint_residual"]) and 0 <= row["maximum_scaled_midpoint_residual"] <= 1.0 for row in rows),
        "svd_equation_residuals": all(below(row["maximum_svd_residual"], 2e-8) for row in rows),
        "fine_norm_conservation": all(below(frame["refinement"][-1]["maximum_norm_drift"], 2e-7) for frame in (fixed, moving)),
        "gauge_error_refines": all(finite(x) and x > 0 for x in gauges) and gauges[0] > 8*gauges[-1],
        "fine_gauge_wavefunction_agreement": below(gauges[-1], 2e-5),
        "physical_triplet_gauge_agreement": below(evidence["fine_triplet_gauge_error"], 2e-6),
        "nonzero_soc_transfer": 1e-3 < fixed["refinement"][-1]["triplet_population"] < 1.0,
        "triplet_agrees_with_grid_reference": all(below(frame["refinement"][-1]["triplet_grid_error"], 5e-4) for frame in (fixed, moving)),
        "wavefunction_agrees_with_grid_reference": all(below(frame["refinement"][-1]["grid_exact_wavefunction_error"], 1e-2) for frame in (fixed, moving)),
        "disabled_soc_blocks_transfer": all(below(value, 1e-12) for value in (off["triplet_population"], off["grid_exact_triplet_population"])),
        "grid_reference_gauge_covariance": below(evidence["grid_exact_gauge_error"], 2e-10),
        "same_runtime_exact_replay": evidence["same_runtime_exact_replay"] is True,
    }.items()}


def _run_soc_frame_v303(moving, verbose=False):
    final_time, step_counts = 0.4, (2, 4, 8)
    progress = (lambda message: print(message, flush=True)) if verbose else None
    key = "moving" if moving else "fixed"
    state, bridge, grid, gauge, config = soc_packet_fixture_v303(moving=moving)
    initial = build_molecular_tdvp_snapshot_v301(state, bridge, grid).psi
    exact = expm_multiply(-1j * final_time * grid.matrix, initial.reshape(-1)).reshape(initial.shape)
    exact_physical = _physical_wavefunction(exact, grid, gauge, moving)
    exact_triplet = _triplet_population(exact_physical, grid)
    frame_endpoints, rows = [], []
    for steps in step_counts:
        if progress is not None:
            progress(f"{key} frame: {steps} steps to t={final_time}")
        trajectory = run_molecular_midpoint_v302(state, bridge, grid, dt=final_time/steps, steps=steps)
        psi = build_molecular_tdvp_snapshot_v301(trajectory[-1].after, bridge, grid).psi
        physical = _physical_wavefunction(psi, grid, gauge, moving)
        frame_endpoints.append(physical)
        triplet = _triplet_population(physical, grid)
        first = trajectory[0].receipt["before_observables"]
        rows.append({
            "dt": final_time/steps, "steps": steps,
            "maximum_norm_drift": max(abs(s.receipt["after_observables"]["norm"] - first["norm"]) for s in trajectory),
            "maximum_energy_drift_hartree": max(abs(s.receipt["after_observables"]["energy_hartree"] - first["energy_hartree"]) for s in trajectory),
            "grid_exact_wavefunction_error": _norm(psi-exact, grid),
            "triplet_population": triplet, "triplet_grid_error": abs(triplet-exact_triplet),
            "maximum_scaled_midpoint_residual": max(s.receipt["scaled_midpoint_residual"] for s in trajectory),
            "maximum_svd_residual": max(s.receipt["midpoint_svd_residual"] for s in trajectory),
        })
        if progress is not None:
            progress(f"{key} result: {json.dumps(rows[-1], sort_keys=True)}")
        if steps == step_counts[0]:
            coarse_trajectory = json.loads(serialize_molecular_trajectory_v302(trajectory))
            if moving:
                if progress is not None:
                    progress("replaying the complete coarse moving-frame trajectory")
                serialized = serialize_molecular_trajectory_v302(trajectory)
                replay_exact = serialize_molecular_trajectory_v302(replay_molecular_trajectory_v302(serialized, bridge, grid)) == serialized
    differences = [_norm(frame_endpoints[i]-frame_endpoints[i+1], grid) for i in range(2)]
    frame = {
        "operator": grid.evidence(), "refinement": rows,
        "successive_wavefunction_differences": differences,
        "observed_time_order": _order(differences),
        "grid_exact_triplet_population": exact_triplet,
    }
    return {"frame": frame, "endpoints": frame_endpoints,
            "exact": exact_physical, "weights": grid.quadrature.weights,
            "config": asdict(config), "trajectory": coarse_trajectory,
            "replay_exact": replay_exact if moving else None}

def _run_soc_disabled_v303(verbose=False):
    final_time = 0.4
    if verbose:
        print("SOC-disabled control: 8 steps", flush=True)
    state, bridge, grid, _, _ = soc_packet_fixture_v303(soc_enabled=False)
    off_initial = build_molecular_tdvp_snapshot_v301(state, bridge, grid).psi
    off_exact = expm_multiply(-1j*final_time*grid.matrix, off_initial.reshape(-1)).reshape(off_initial.shape)
    off_trajectory = run_molecular_midpoint_v302(state, bridge, grid, dt=final_time/8, steps=8)
    off_psi = build_molecular_tdvp_snapshot_v301(off_trajectory[-1].after, bridge, grid).psi
    return {"triplet_population": _triplet_population(off_psi, grid),
            "grid_exact_triplet_population": _triplet_population(off_exact, grid)}


def run_v303_soc_campaign(*, workers=1, verbose=False):
    """Propagate to T=0.4 a.u.; independent cases may use separate processes.

    Each process owns its provider, random seed and complete trajectory/replay.
    Results are collected in a fixed order; process scheduling enters no receipt.
    """
    if type(workers) is not int or workers not in (1, 2, 3):
        raise ValueError("workers must be 1, 2 or 3")
    if workers == 1:
        fixed = _run_soc_frame_v303(False, verbose)
        moving = _run_soc_frame_v303(True, verbose)
        off = _run_soc_disabled_v303(verbose)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            fixed_job = executor.submit(_run_soc_frame_v303, False, verbose)
            moving_job = executor.submit(_run_soc_frame_v303, True, verbose)
            off_job = executor.submit(_run_soc_disabled_v303, verbose)
            fixed, moving, off = fixed_job.result(), moving_job.result(), off_job.result()
    weights = fixed["weights"]
    if not np.array_equal(weights, moving["weights"]):
        raise ValueError("frame quadratures disagree")
    frames = {"fixed": fixed["frame"], "moving": moving["frame"]}
    evidence = {
        "schema": SCHEMA_V303, "source_kind": "analytic complete singlet-triplet model",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "model_config": fixed["config"], "final_time_au": 0.4,
        "frames": frames,
        "gauge_wavefunction_errors": [_weighted_norm(a-b, weights) for a,b in zip(fixed["endpoints"], moving["endpoints"])],
        "fine_triplet_gauge_error": abs(frames["fixed"]["refinement"][-1]["triplet_population"]-frames["moving"]["refinement"][-1]["triplet_population"]),
        "grid_exact_gauge_error": _weighted_norm(fixed["exact"]-moving["exact"], weights),
        "soc_disabled": off,
        "same_runtime_exact_replay": moving["replay_exact"],
        "coarse_trajectories": {"fixed": fixed["trajectory"], "moving": moving["trajectory"]},
        "live_molecular_soc_validated": False, "dataset_certified": False,
        "scope": "One packet; analytic 1D model; fixed finite grid. Temporal and gauge convergence do not establish continuum or basis completeness.",
    }
    evidence["gates"] = evaluate_soc_gates_v303(evidence)
    evidence["passed"] = all(evidence["gates"].values())
    return evidence
