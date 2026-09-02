from gaussian_dynamics.v252_benchmark import run_v0252_release_benchmark


def test_v0252_release_campaign_passes_630_native_boolean_gates():
    output = run_v0252_release_benchmark()
    acceptance = output["acceptance"]

    assert output["release"] == "v0.25.2"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 535
    assert acceptance["validation_gate_count"] == 70
    assert acceptance["core_gate_count"] == 25
    assert acceptance["new_gate_count"] == 95
    assert acceptance["total_gate_count"] == 630
    assert len(acceptance["checks"]) == 630
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
