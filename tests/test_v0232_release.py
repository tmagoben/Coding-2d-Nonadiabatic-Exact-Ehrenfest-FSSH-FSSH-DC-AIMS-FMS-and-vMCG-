from gaussian_dynamics import run_v0232_release_benchmark


def test_v0232_release_campaign_passes_168_native_boolean_gates():
    output = run_v0232_release_benchmark()

    acceptance = output["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 123
    assert acceptance["runtime_gate_count"] == 28
    assert acceptance["new_gate_count"] == 45
    assert acceptance["total_gate_count"] == 168
    assert len(acceptance["checks"]) == 168
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = output["claims"]
    assert claims["real_PySCF_spin_free_runtime_validated"] is True
    assert claims["real_PySCF_SA_CASSCF_gradients_validated"] is True
    assert claims["real_PySCF_NAC_overlap_consistency_validated"] is True
    assert claims["finite_manifold_overlap_contract_validated"] is True
    assert claims["trust_anchored_runtime_admission_validated"] is True
    assert claims["physical_analytic_SOC_inherited"] is True
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["ab_initio_SOC_validated"] is False
    assert claims["live_PySCF_SOC_runtime_validated"] is False
