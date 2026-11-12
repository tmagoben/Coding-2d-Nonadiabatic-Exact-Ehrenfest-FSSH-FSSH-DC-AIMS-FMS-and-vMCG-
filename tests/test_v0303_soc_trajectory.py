"""SOC campaign contracts; CI separately recomputes the complete trajectories."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from gaussian_dynamics.complex_gauge_v21 import PhaseMixingGaugeV21, random_unitary_v21
from gaussian_dynamics.curved_basis_adaptation_v292 import CurvedQuadratureV292
from gaussian_dynamics.v303_benchmark import (
    SCHEMA_V303, _physical_wavefunction, _triplet_population,
    evaluate_soc_gates_v303,
)


@pytest.fixture
def evidence():
    return json.loads((Path(__file__).resolve().parents[1] / "results/v0303_soc_trajectory_evidence.json").read_text())


def test_committed_evidence_is_consistent_and_scope_is_explicit(evidence):
    assert evidence["schema"] == SCHEMA_V303
    assert evidence["gates"] == evaluate_soc_gates_v303(evidence)
    assert evidence["passed"] and all(evidence["gates"].values())
    assert not evidence["live_molecular_soc_validated"]
    assert not evidence["dataset_certified"]
    assert evidence["final_time_au"] == 0.4
    for frame in evidence["frames"].values():
        assert [row["steps"] for row in frame["refinement"]] == [2, 4, 8]
        assert all(row["dt"]*row["steps"] == 0.4 for row in frame["refinement"])
        # Report the remaining projection error even after timestep convergence.
        assert frame["refinement"][-1]["grid_exact_wavefunction_error"] > frame["successive_wavefunction_differences"][-1]


@pytest.mark.parametrize("frame", ["fixed", "moving"])
@pytest.mark.parametrize("field,bad,gate", [
    ("maximum_scaled_midpoint_residual", 1.1, "accepted_midpoint_residuals"),
    ("maximum_svd_residual", 1e-5, "svd_equation_residuals"),
    ("maximum_norm_drift", 1e-4, "fine_norm_conservation"),
    ("triplet_grid_error", 0.01, "triplet_agrees_with_grid_reference"),
    ("grid_exact_wavefunction_error", 0.1, "wavefunction_agrees_with_grid_reference"),
])
def test_failed_numerics_cannot_be_overridden_by_pass_flags(evidence, frame, field, bad, gate):
    evidence["frames"][frame]["refinement"][-1][field] = bad
    evidence["gates"] = {key: True for key in evidence["gates"]}
    assert not evaluate_soc_gates_v303(evidence)[gate]


@pytest.mark.parametrize("bad", [0.0, 1.0, np.nan, np.inf])
def test_false_time_convergence_is_rejected(evidence, bad):
    evidence["frames"]["moving"]["observed_time_order"] = bad
    assert not evaluate_soc_gates_v303(evidence)["moving_frame_second_order"]


def test_claimed_order_must_match_the_wavefunction_differences(evidence):
    evidence["frames"]["fixed"]["successive_wavefunction_differences"] = [1e-4, 1e-4]
    assert not evaluate_soc_gates_v303(evidence)["fixed_frame_second_order"]


def test_missing_soc_signal_and_disabled_control_contamination_fail(evidence):
    no_signal = deepcopy(evidence)
    no_signal["frames"]["fixed"]["refinement"][-1]["triplet_population"] = evidence["soc_disabled"]["triplet_population"]
    assert not evaluate_soc_gates_v303(no_signal)["nonzero_soc_transfer"]
    evidence["soc_disabled"]["triplet_population"] = 1e-3
    assert not evaluate_soc_gates_v303(evidence)["disabled_soc_blocks_transfer"]


def test_unresolved_gauge_error_and_failed_replay_are_rejected(evidence):
    evidence["gauge_wavefunction_errors"] = [1e-5, 1e-5, 1e-5]
    evidence["same_runtime_exact_replay"] = False
    gates = evaluate_soc_gates_v303(evidence)
    assert not gates["gauge_error_refines"]
    assert not gates["same_runtime_exact_replay"]


def test_spin_populations_require_transforming_out_of_the_mixed_gauge():
    from types import SimpleNamespace
    points = np.array([[-0.4], [0.0], [0.4]])
    grid = SimpleNamespace(quadrature=CurvedQuadratureV292(points, np.full(3, 0.4)))
    gauge = PhaseMixingGaugeV21(random_unitary_v21(4, seed=303), np.array([[.3],[-.2],[.1],[.25]]), np.zeros(4))
    fixed = np.array([[1.0, 0, 0, 0]]*3, dtype=complex)
    moving = np.array([gauge.matrix(q).conj().T @ psi for q, psi in zip(points, fixed)])
    observed = _physical_wavefunction(moving, grid, gauge, True)
    np.testing.assert_allclose(observed, fixed, atol=1e-14)
    assert _triplet_population(observed, grid) < 1e-28
    # Negative control: taking mixed-gauge component magnitudes as triplets
    # invents a large transition even for a physically pure singlet.
    assert _triplet_population(moving, grid) > 0.1
