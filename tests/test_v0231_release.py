from gaussian_dynamics import run_v0231_release_benchmark


def test_v0231_release_campaign_passes_123_native_boolean_gates():
    output = run_v0231_release_benchmark()

    assert output["acceptance"]["passed"]
    assert output["acceptance"]["inherited_gate_count"] == 93
    assert output["acceptance"]["new_gate_count"] == 30
    assert output["acceptance"]["total_gate_count"] == 123
    assert len(output["acceptance"]["checks"]) == 123
    assert all(type(value) is bool for value in output["acceptance"]["checks"].values())
    assert all(output["acceptance"]["checks"].values())
    assert output["claims"]["raw_evidence_admission_protocol_validated"] is True
    assert output["claims"]["external_molecular_SOC_snapshot_admitted"] is False
    assert output["claims"]["live_molecular_SOC_backend_admitted"] is False
    assert output["claims"]["ab_initio_SOC_validated"] is False
