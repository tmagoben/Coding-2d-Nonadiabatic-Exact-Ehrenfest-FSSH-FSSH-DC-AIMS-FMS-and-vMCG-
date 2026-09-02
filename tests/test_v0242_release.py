from gaussian_dynamics.v242_benchmark import run_v0242_release_benchmark


def test_v0242_release_campaign_passes_400_native_boolean_gates():
    output = run_v0242_release_benchmark()
    acceptance = output["acceptance"]

    assert output["release"] == "v0.24.2"
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 315
    assert acceptance["runtime_gate_count"] == 60
    assert acceptance["core_gate_count"] == 25
    assert acceptance["new_gate_count"] == 85
    assert acceptance["total_gate_count"] == 400
    assert len(acceptance["checks"]) == 400
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())
