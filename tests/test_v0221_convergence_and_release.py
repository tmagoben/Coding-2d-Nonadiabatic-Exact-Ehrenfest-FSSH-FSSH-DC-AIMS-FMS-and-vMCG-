from gaussian_dynamics import run_v0221_release_benchmark
from gaussian_dynamics.v221_benchmark import (
    _basis_convergence_v221,
    _sparse_threshold_convergence_v221,
)


def test_soc_gaussian_basis_population_ladder_narrows():
    result = _basis_convergence_v221()

    assert [row["basis_size"] for row in result["rows"]] == [1, 3, 5]
    assert result["coarse_population_difference"] > result["fine_population_difference"]
    assert result["fine_population_difference"] < 1.0e-8
    assert result["narrowing_ratio"] < 0.05


def test_soc_sparse_threshold_ladder_converges_to_dense():
    result = _sparse_threshold_convergence_v221()
    errors = [row["coefficient_error"] for row in result["rows"]]
    edges = [row["active_edges"] for row in result["rows"]]

    assert all(later <= earlier + 1.0e-14 for earlier, later in zip(errors, errors[1:]))
    assert all(later >= earlier for earlier, later in zip(edges, edges[1:]))
    assert errors[-1] < 1.0e-12


def test_v0221_release_campaign_passes_67_native_boolean_gates():
    output = run_v0221_release_benchmark()

    assert output["acceptance"]["passed"]
    assert output["acceptance"]["inherited_gate_count"] == 53
    assert output["acceptance"]["new_gate_count"] == 14
    assert output["acceptance"]["total_gate_count"] == 67
    assert len(output["acceptance"]["checks"]) == 67
    assert all(type(value) is bool for value in output["acceptance"]["checks"].values())
    assert all(output["acceptance"]["checks"].values())
    assert output["inherited_v0220"]["acceptance"]["passed"]
    assert output["soc"]["analytic_models_only"] is True
    assert output["soc"]["ab_initio_SOC_validated"] is False
    assert output["soc"]["molecular_SOC_backend_admitted"] is False

