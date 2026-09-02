from dataclasses import replace
import copy

import numpy as np
import pytest

from gaussian_dynamics import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    SOCOperatorComponentsV220,
    SOCSymmetryContractV221,
    SpinorGridSettingsV220,
    audit_physical_soc_provider_v220,
    audit_soc_symmetry_contract_v221,
    initial_gaussian_spinor_v220,
    phase_aligned_spinor_grid_error_v220,
    projector_population_v220,
    run_spinor_exact_grid_v220,
    spinor_split_operator_step_v220,
    transform_projector_v220,
)
from gaussian_dynamics.v221_benchmark import (
    _CancelledComponentDerivativeV221,
    _GenericThreeStateTwoCoordinateProviderV221,
    _NoConfigProviderV221,
    _grid_hardening_v221,
    _symmetry_hardening_v221,
)


def test_component_container_is_dimension_neutral_and_canonicalizes_inputs():
    components = SOCOperatorComponentsV220(
        [0.1, -0.2],
        np.eye(3).tolist(),
        np.zeros((2, 3, 3)).tolist(),
        np.zeros((3, 3)).tolist(),
        np.zeros((2, 3, 3)).tolist(),
    ).validate()

    assert components.q.shape == (2,)
    assert components.H.shape == (3, 3)
    assert components.K.shape == (2, 3, 3)
    assert components.H.dtype == complex


def test_three_state_two_coordinate_no_config_provider_passes_full_audit():
    provider = _GenericThreeStateTwoCoordinateProviderV221()
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([0.13, -0.21]), fermionic=False
    )

    assert not hasattr(provider, "config")
    assert report.passed
    assert report.maximum_spin_free_component_derivative_error < 2.0e-9
    assert report.maximum_soc_component_derivative_error < 2.0e-9
    assert len(report.component_derivative_rows) == 6


def test_cancelled_spin_free_and_soc_derivative_errors_are_rejected():
    provider = _CancelledComponentDerivativeV221(
        AnalyticSingletTripletSOCProviderV220()
    )
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([0.17]), fermionic=False
    )

    assert not report.passed
    assert report.checks["K_decomposition"]
    assert report.checks["cross_geometry_differentials"]
    assert report.checks["SOC_force_derivative"]
    assert not report.checks["spin_free_component_derivatives"]
    assert not report.checks["SOC_component_derivatives"]


def test_provider_config_is_not_part_of_the_soc_audit_contract():
    provider = _NoConfigProviderV221(AnalyticDoubletSOCProviderV220())
    report = audit_physical_soc_provider_v220(
        provider, np.asarray([-0.11]), fermionic=True
    )

    assert not hasattr(provider, "config")
    assert report.passed


def test_mixed_electron_parity_is_rejected():
    report = _symmetry_hardening_v221()["mixed_parity"]

    assert not report["passed"]
    assert not report["checks"]["single_electron_parity"]


def test_time_reversal_unitarity_is_independent_of_its_square():
    report = _symmetry_hardening_v221()["nonunitary_time_reversal"]

    assert report["time_reversal_square_residual"] == 0.0
    assert report["time_reversal_unitarity_residual"] > 0.5
    assert not report["checks"]["time_reversal_unitarity"]
    assert not report["passed"]


def test_numerical_time_reversal_and_projectors_are_provenance_identity():
    provider = AnalyticSingletTripletSOCProviderV220()
    valid = audit_soc_symmetry_contract_v221(
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
        provenance=provider.provenance,
        fermionic=False,
    )
    parameters = copy.deepcopy(provider.provenance.parameters)
    parameters["physical_projectors"] = {
        "singlet": provider.projectors["triplet"].tolist(),
        "triplet": provider.projectors["singlet"].tolist(),
    }
    changed = replace(provider.provenance, parameters=parameters)
    mismatch = audit_soc_symmetry_contract_v221(
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
        provenance=changed,
        fermionic=False,
    )

    assert valid.passed
    assert not mismatch.passed
    assert not mismatch.checks["symmetry_provenance_identity"]


def test_nonunitary_gauge_and_invalid_projector_fail_closed():
    bad_gauge = np.asarray([[1.0, 1.0], [0.0, 1.0]], dtype=complex)
    with pytest.raises(ValueError, match="unitary"):
        transform_projector_v220(np.diag([1.0, 0.0]), bad_gauge)
    with pytest.raises(ValueError, match="Hermitian projector"):
        projector_population_v220(
            np.asarray([1.0, 0.0]), np.asarray([[1.0, 1.0], [0.0, 0.0]])
        )


def test_exact_grid_always_records_the_actual_final_step():
    report = _grid_hardening_v221()

    assert report["recorded_times"] == [0.0, 0.02, 0.04, 0.05]
    assert report["recorded_final_time"] == report["requested_final_time"]


def test_exact_grid_uses_emitted_mass_without_provider_config():
    x = np.linspace(-4.0, 4.0, 64, endpoint=False)
    provider = _NoConfigProviderV221(AnalyticSingletTripletSOCProviderV220())
    psi0 = initial_gaussian_spinor_v220(
        x, np.asarray([1.0, 0.0, 0.0, 0.0])
    )
    output = run_spinor_exact_grid_v220(
        provider,
        x,
        psi0,
        settings=SpinorGridSettingsV220(dt=0.01, steps=2, store_every=2),
    )

    assert not hasattr(provider, "config")
    assert output["constant_mass_certified"]
    assert output["fixed_frame_certified"]


def test_exact_grid_rejects_a_moving_electronic_frame():
    report = _grid_hardening_v221()
    assert report["moving_frame_rejected"]


def test_precomputed_grid_run_matches_public_step_operator():
    x = np.linspace(-5.0, 5.0, 128, endpoint=False)
    dx = x[1] - x[0]
    provider = AnalyticSingletTripletSOCProviderV220()
    psi0 = initial_gaussian_spinor_v220(
        x, np.asarray([1.0, 0.0, 0.0, 0.0])
    )
    output = run_spinor_exact_grid_v220(
        provider,
        x,
        psi0,
        settings=SpinorGridSettingsV220(dt=0.02, steps=3, store_every=3),
    )
    potential = np.asarray(
        [provider.evaluate_snapshot(np.asarray([coordinate])).point.H for coordinate in x]
    )
    manual = psi0.copy()
    for _ in range(3):
        manual = spinor_split_operator_step_v220(
            manual, dx, 0.02, 900.0, potential
        )

    assert phase_aligned_spinor_grid_error_v220(
        output["psi"][-1], manual, dx
    ) < 1.0e-14


def test_grid_error_rejects_nonpositive_spacing():
    state = np.ones((2, 8), dtype=complex)
    with pytest.raises(ValueError, match="spacing"):
        phase_aligned_spinor_grid_error_v220(state, state, 0.0)
