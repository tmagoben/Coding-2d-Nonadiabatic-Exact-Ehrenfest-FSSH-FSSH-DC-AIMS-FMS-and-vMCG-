from gaussian_dynamics.v251_benchmark import run_v0251_release_benchmark


def test_v0251_release_campaign_passes_535_native_boolean_gates():
    output = run_v0251_release_benchmark()
    acceptance = output["acceptance"]

    assert output["release"] == "v0.25.1"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 460
    assert acceptance["validation_gate_count"] == 55
    assert acceptance["core_gate_count"] == 20
    assert acceptance["new_gate_count"] == 75
    assert acceptance["total_gate_count"] == 535
    assert len(acceptance["checks"]) == 535
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
