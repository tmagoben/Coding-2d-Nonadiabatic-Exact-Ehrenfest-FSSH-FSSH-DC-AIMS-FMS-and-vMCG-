from gaussian_dynamics import run_v0233_release_benchmark


def test_v0233_release_campaign_passes_208_native_boolean_gates():
    output = run_v0233_release_benchmark()

    acceptance = output["acceptance"]
    assert acceptance["passed"] is True
    assert acceptance["inherited_gate_count"] == 168
    assert acceptance["new_gate_count"] == 40
    assert acceptance["total_gate_count"] == 208
    assert len(acceptance["checks"]) == 208
    assert all(type(value) is bool for value in acceptance["checks"].values())
    assert all(acceptance["checks"].values())

    claims = output["claims"]
    assert claims["finite_manifold_unitary_transport_validated"] is True
    assert claims["replay_format_two_and_migration_validated"] is True
    assert claims["legacy_NAC_data_quarantine_validated"] is True
    assert claims["complete_even_and_odd_manifold_transport_validated"] is True
    assert claims["molecular_SOC_matrix_convention_frozen"] is True
    assert claims["runtime_identity_and_compatibility_profiles_separated"] is True
    assert claims["real_PySCF_spin_free_runtime_inherited"] is True
    assert claims["external_molecular_SOC_snapshot_admitted"] is False
    assert claims["live_molecular_SOC_backend_admitted"] is False
    assert claims["ab_initio_SOC_validated"] is False
    assert claims["live_PySCF_SOC_runtime_validated"] is False
