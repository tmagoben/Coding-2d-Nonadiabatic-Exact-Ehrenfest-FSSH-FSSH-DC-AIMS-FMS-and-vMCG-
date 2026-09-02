from gaussian_dynamics.v250_benchmark import run_v0250_release_benchmark


def test_v0250_release_campaign_passes_460_native_boolean_gates():
    output = run_v0250_release_benchmark()
    acceptance = output["acceptance"]

    assert output["release"] == "v0.25.0"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 400
    assert acceptance["validation_gate_count"] == 45
    assert acceptance["core_gate_count"] == 15
    assert acceptance["new_gate_count"] == 60
    assert acceptance["total_gate_count"] == 460
    assert len(acceptance["checks"]) == 460
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
