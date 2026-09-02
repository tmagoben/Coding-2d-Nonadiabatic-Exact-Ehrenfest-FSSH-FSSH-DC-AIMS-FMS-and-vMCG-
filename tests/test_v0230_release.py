from gaussian_dynamics import run_v0230_release_benchmark


def test_v0230_release_campaign_passes_93_native_boolean_gates():
    output = run_v0230_release_benchmark()

    assert output["acceptance"]["passed"]
    assert output["acceptance"]["inherited_gate_count"] == 67
    assert output["acceptance"]["new_gate_count"] == 26
    assert output["acceptance"]["total_gate_count"] == 93
    assert len(output["acceptance"]["checks"]) == 93
    assert all(type(value) is bool for value in output["acceptance"]["checks"].values())
    assert all(output["acceptance"]["checks"].values())
    assert output["claims"]["molecular_SOC_protocol_validated"] is True
    assert output["claims"]["real_molecular_SOC_backend_admitted"] is False
    assert output["claims"]["ab_initio_SOC_validated"] is False
