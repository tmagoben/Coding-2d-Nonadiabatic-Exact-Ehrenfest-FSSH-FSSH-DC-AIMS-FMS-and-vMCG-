from gaussian_dynamics import run_v0241_release_benchmark


def test_v0241_release_campaign_passes_315_native_boolean_gates():
    output = run_v0241_release_benchmark()
    acceptance = output["acceptance"]

    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 256
    assert acceptance["runtime_gate_count"] == 39
    assert acceptance["core_gate_count"] == 20
    assert acceptance["new_gate_count"] == 59
    assert acceptance["total_gate_count"] == 315
    assert len(acceptance["checks"]) == 315
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = output["claims"]
    assert claims["real_PySCF_BP_SOMF_execution_validated"] is True
    assert claims["direct_molecular_SOC_elements_returned"] is True
    assert claims["doublet_and_Kramers_sector_validated"] is True
    assert claims["mixed_multiplicity_spin_algebra_validated"] is True
    assert claims["static_molecular_SOC_tier_validated"] is True
    assert claims["trajectory_ready_molecular_SOC_validated"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["physical_SOC_derivatives_validated"] is False
    assert claims["cross_geometry_SOC_tracking_validated"] is False
    assert claims["ab_initio_SOC_accuracy_validated"] is False
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["native_openmolcas_numeric_crosscheck_implemented"] is False
    assert claims["Prism_runtime_dependency_required"] is False
