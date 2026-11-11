from gaussian_dynamics.v270_benchmark import (
    V270AcceptanceThresholds,
    run_v0270_release_benchmark,
)


def test_v0270_cumulative_release_campaign_checks_all_gates(assert_runtime_aware_campaign):
    result = run_v0270_release_benchmark()
    acceptance = result["acceptance"]
    thresholds = V270AcceptanceThresholds()
    assert result["release"] == "v0.27.0"
    assert_runtime_aware_campaign(result)
    assert acceptance["inherited_gate_count"] == thresholds.expected_inherited_gates
    assert acceptance["validation_gate_count"] == thresholds.expected_validation_gates
    assert acceptance["core_gate_count"] == thresholds.expected_core_gates
    assert acceptance["new_gate_count"] == thresholds.expected_new_gates
    assert acceptance["total_gate_count"] == thresholds.expected_total_gates
    assert len(acceptance["checks"]) == 960
