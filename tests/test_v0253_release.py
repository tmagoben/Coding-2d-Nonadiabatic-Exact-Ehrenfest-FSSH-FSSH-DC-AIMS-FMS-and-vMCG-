from gaussian_dynamics.v253_benchmark import (
    V253AcceptanceThresholds,
    run_v0253_release_benchmark,
)


def test_v0253_cumulative_release_campaign_passes_all_gates():
    result = run_v0253_release_benchmark()
    acceptance = result["acceptance"]
    thresholds = V253AcceptanceThresholds()
    assert result["release"] == "v0.25.3"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == thresholds.expected_inherited_gates
    assert acceptance["validation_gate_count"] == thresholds.expected_validation_gates
    assert acceptance["core_gate_count"] == thresholds.expected_core_gates
    assert acceptance["new_gate_count"] == thresholds.expected_new_gates
    assert acceptance["total_gate_count"] == thresholds.expected_total_gates
    assert len(acceptance["checks"]) == 715
    assert all(acceptance["checks"].values())
