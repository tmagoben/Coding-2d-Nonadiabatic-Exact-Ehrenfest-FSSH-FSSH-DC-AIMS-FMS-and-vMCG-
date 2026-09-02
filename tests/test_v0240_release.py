from gaussian_dynamics import run_v0240_release_benchmark


def test_v0240_release_campaign_passes_256_native_boolean_gates():
    output = run_v0240_release_benchmark()
    acceptance = output["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 208
    assert acceptance["new_gate_count"] == 48
    assert acceptance["total_gate_count"] == 256
    assert len(acceptance["checks"]) == 256
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = output["claims"]
    assert claims["openmolcas_rassi_so_protocol_frozen"] is True
    assert claims["strict_bundle_artifact_parser_validated"] is True
    assert claims["transported_cartesian_soc_derivative_protocol_validated"] is True
    assert claims["protocol_fixture_validated"] is True
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["ab_initio_SOC_validated"] is False
    assert claims["openmolcas_runtime_executed"] is False
    assert claims["native_openmolcas_numeric_crosscheck_implemented"] is False
