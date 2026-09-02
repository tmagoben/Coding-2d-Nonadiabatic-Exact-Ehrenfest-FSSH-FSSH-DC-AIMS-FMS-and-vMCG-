from gaussian_dynamics.v260_benchmark import (
    V260AcceptanceThresholds,
    run_v0260_release_benchmark,
)


def test_v0260_cumulative_release_campaign_passes_all_gates():
    result = run_v0260_release_benchmark()
    acceptance = result["acceptance"]
    thresholds = V260AcceptanceThresholds()
    assert result["release"] == "v0.26.0"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == thresholds.expected_inherited_gates
    assert acceptance["validation_gate_count"] == thresholds.expected_validation_gates
    assert acceptance["core_gate_count"] == thresholds.expected_core_gates
    assert acceptance["new_gate_count"] == thresholds.expected_new_gates
    assert acceptance["total_gate_count"] == thresholds.expected_total_gates
    assert len(acceptance["checks"]) == 825
    assert all(acceptance["checks"].values())
