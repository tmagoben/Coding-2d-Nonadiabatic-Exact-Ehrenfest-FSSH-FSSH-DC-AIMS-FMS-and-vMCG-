from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
)
from gaussian_dynamics.complex_gauge_v21 import (
    GaugeTransformedOperatorProviderV21,
    PhaseMixingGaugeV21,
    random_unitary_v21,
)
from gaussian_dynamics.multigaussian_tdvp_v251 import (
    MULTIGAUSSIAN_TDVP_ANSATZ_V251,
    VARIATIONAL_INTEGRATOR_V251,
    VARIATIONAL_METRIC_SOLVER_V251,
    FrozenGaussianSpinorStateV251,
    QuadraticSpinHamiltonianV251,
    VariationalMetricSettingsV251,
    build_frozen_gaussian_spinor_matrices_v251,
    build_variational_metric_system_v251,
    implicit_midpoint_tdvp_step_v251,
    pack_variational_parameters_v251,
    quadratic_spin_hamiltonian_from_provider_v251,
    reverse_frozen_width_multigaussian_tdvp_v251,
    run_frozen_width_multigaussian_tdvp_v251,
    solve_variational_metric_v251,
    state_from_variational_parameters_v251,
)


def _coefficients():
    return np.asarray(
        [
            [0.65 + 0.10j, 0.15 - 0.20j, 0.25 + 0.08j, -0.05j],
            [0.18 - 0.04j, -0.11 + 0.09j, 0.22 - 0.06j, 0.07 + 0.03j],
        ],
        dtype=complex,
    )


def _state():
    return FrozenGaussianSpinorStateV251(
        q=np.asarray([-0.65, 0.75]),
        p=np.asarray([5.0, -3.0]),
        widths=np.asarray([2.6, 2.1]),
        coefficients=_coefficients(),
    ).normalized()


def _model():
    return quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220()
    )


def _maximum_state_error(left, right):
    return max(
        float(np.max(np.abs(left.q - right.q))),
        float(np.max(np.abs(left.p - right.p))),
        float(np.max(np.abs(left.widths - right.widths))),
        float(np.max(np.abs(left.coefficients - right.coefficients))),
        abs(float(left.time_au) - float(right.time_au)),
    )


def test_scope_freezes_metric_midpoint_and_keeps_later_features_closed():
    settings = VariationalMetricSettingsV251().validate()
    assert "multi" in MULTIGAUSSIAN_TDVP_ANSATZ_V251
    assert settings.integrator == VARIATIONAL_INTEGRATOR_V251
    assert settings.metric_solver == VARIATIONAL_METRIC_SOLVER_V251
    assert settings.allow_compatible_rank_deficiency is True

    for field_name in (
        "adaptive_gaussian_widths",
        "spawning",
        "pruning",
        "coordinate_dependent_electronic_frame",
        "multidimensional_nuclear_motion",
        "real_molecular_soc_provider",
    ):
        with pytest.raises(ValueError, match="does not admit"):
            replace(settings, **{field_name: True}).validate()
    with pytest.raises(ValueError, match="SVD metric solver is frozen"):
        replace(settings, metric_solver="normal equations").validate()
    with pytest.raises(ValueError, match="implicit integrator is frozen"):
        replace(settings, integrator="velocity Verlet").validate()


@pytest.mark.parametrize(
    "provider",
    [AnalyticSingletTripletSOCProviderV220(), AnalyticDoubletSOCProviderV220()],
)
def test_quadratic_intake_reproduces_complete_even_and_odd_soc_providers(provider):
    model = quadratic_spin_hamiltonian_from_provider_v251(provider)
    assert model.nstate == 4
    assert model.physical_soc is True
    assert model.complete_spin_manifold is True
    assert model.source["fixed_frame_verified"] is True
    assert model.source["model_space"]["complete_multiplets"] is True
    for x in (-1.71, -0.24, 0.63, 1.48):
        point = provider.evaluate_snapshot(np.asarray([x])).point
        assert np.allclose(model.hamiltonian(x), point.H, atol=3.0e-14)
        assert np.allclose(
            model.derivative(x),
            point.hamiltonian_derivative_operator_q[0],
            atol=3.0e-14,
        )


def test_quadratic_intake_rejects_moving_frames_and_static_receipts():
    base = AnalyticDoubletSOCProviderV220()
    gauge = PhaseMixingGaugeV21(
        random_unitary_v21(4, 25102),
        np.asarray([[0.17], [-0.11], [0.23], [-0.07]]),
        np.asarray([0.20, -0.30, 0.10, 0.40]),
    )
    with pytest.raises(TypeError, match="explicit operator provenance"):
        quadratic_spin_hamiltonian_from_provider_v251(
            GaugeTransformedOperatorProviderV21(base, gauge)
        )
    static = SimpleNamespace(evaluate_snapshot=lambda q: SimpleNamespace(matrices=1))
    with pytest.raises(TypeError, match="explicit operator provenance"):
        quadratic_spin_hamiltonian_from_provider_v251(static)


def test_state_normalization_pack_roundtrip_and_validation_are_explicit():
    raw = FrozenGaussianSpinorStateV251(
        [-0.4, 0.6], [1.0, -2.0], [2.4, 2.0], 2.0 * _coefficients()
    )
    with pytest.raises(ValueError, match="unit generalized norm"):
        raw.validate()
    state = raw.normalized()
    assert state.generalized_norm == pytest.approx(1.0, abs=3.0e-16)
    packed = pack_variational_parameters_v251(state)
    restored = state_from_variational_parameters_v251(
        packed, widths=state.widths, nstate=state.nstate, time_au=state.time_au
    )
    assert _maximum_state_error(restored, state) == 0.0
    with pytest.raises(ValueError, match="Gaussian permutation"):
        state.permuted(np.asarray([0.0, 1.0]))
    with pytest.raises(ValueError, match="widths must be positive"):
        FrozenGaussianSpinorStateV251(
            [-0.4], [1.0], [-2.0], [[1.0, 0.0, 0.0, 0.0]]
        ).validate(require_normalized=False)


def test_exact_gaussian_spinor_matrices_match_independent_grid_quadrature():
    model = QuadraticSpinHamiltonianV251(
        mass_au=7.0,
        H0=np.asarray([[0.04, 0.01 + 0.02j], [0.01 - 0.02j, 0.07]]),
        H1=np.asarray([[0.02, -0.008j], [0.008j, -0.01]]),
        H2=np.asarray([[0.006, 0.002], [0.002, 0.009]]),
    ).validate()
    state = FrozenGaussianSpinorStateV251(
        q=[-0.45, 0.62],
        p=[0.38, -0.27],
        widths=[1.35, 1.75],
        coefficients=[[0.8, 0.1j], [0.2 - 0.1j, 0.3]],
    ).normalized()
    analytic_overlap, analytic_hamiltonian = (
        build_frozen_gaussian_spinor_matrices_v251(state, model)
    )

    x = np.linspace(-9.0, 9.0, 60001)
    gaussians = []
    h_on_gaussians = []
    for q, p, width in zip(state.q, state.p, state.widths):
        y = x - q
        gaussian = (width / np.pi) ** 0.25 * np.exp(
            -0.5 * width * y**2 + 1.0j * p * y
        )
        second_factor = (-width * y + 1.0j * p) ** 2 - width
        kinetic = -0.5 * second_factor * gaussian / model.mass_au
        gaussians.append(gaussian)
        h_on_gaussians.append(kinetic)
    numerical_overlap = np.zeros_like(analytic_overlap)
    numerical_hamiltonian = np.zeros_like(analytic_hamiltonian)
    for i in range(state.ngaussian):
        si = slice(i * state.nstate, (i + 1) * state.nstate)
        for j in range(state.ngaussian):
            sj = slice(j * state.nstate, (j + 1) * state.nstate)
            nuclear_overlap = np.trapezoid(
                np.conj(gaussians[i]) * gaussians[j], x
            )
            numerical_overlap[si, sj] = nuclear_overlap * np.eye(state.nstate)
            kinetic = np.trapezoid(
                np.conj(gaussians[i]) * h_on_gaussians[j], x
            )
            potential = np.zeros((state.nstate, state.nstate), dtype=complex)
            weights = np.conj(gaussians[i]) * gaussians[j]
            for a in range(state.nstate):
                for b in range(state.nstate):
                    values = (
                        model.H0[a, b]
                        + x * model.H1[a, b]
                        + x**2 * model.H2[a, b]
                    )
                    potential[a, b] = np.trapezoid(weights * values, x)
            numerical_hamiltonian[si, sj] = (
                kinetic * np.eye(state.nstate) + potential
            )
    assert np.allclose(analytic_overlap, numerical_overlap, atol=2.0e-11)
    assert np.allclose(analytic_hamiltonian, numerical_hamiltonian, atol=2.0e-11)


def test_analytic_tdvp_metric_and_rhs_match_independent_grid_tangents():
    model = QuadraticSpinHamiltonianV251(
        mass_au=8.0,
        H0=np.asarray([[0.03, 0.012 + 0.008j], [0.012 - 0.008j, 0.06]]),
        H1=np.asarray([[0.015, -0.004j], [0.004j, -0.009]]),
        H2=np.asarray([[0.005, 0.0015], [0.0015, 0.008]]),
    ).validate()
    state = FrozenGaussianSpinorStateV251(
        q=[-0.42, 0.58],
        p=[0.31, -0.22],
        widths=[1.4, 1.8],
        coefficients=[[0.72 + 0.05j, 0.08 - 0.12j], [0.19j, 0.31 - 0.04j]],
    ).normalized()
    analytic = build_variational_metric_system_v251(state, model)
    parameters = pack_variational_parameters_v251(state)
    x = np.linspace(-9.0, 9.0, 40001)

    def wavefunction(theta):
        local = state_from_variational_parameters_v251(
            theta, widths=state.widths, nstate=state.nstate, time_au=0.0
        )
        value = np.zeros((len(x), state.nstate), dtype=complex)
        for q, p, width, coefficients in zip(
            local.q, local.p, local.widths, local.coefficients
        ):
            y = x - q
            gaussian = (width / np.pi) ** 0.25 * np.exp(
                -0.5 * width * y**2 + 1.0j * p * y
            )
            value += gaussian[:, None] * coefficients[None, :]
        return value

    psi = wavefunction(parameters)
    hpsi = np.zeros_like(psi)
    for q, p, width, coefficients in zip(
        state.q, state.p, state.widths, state.coefficients
    ):
        y = x - q
        gaussian = (width / np.pi) ** 0.25 * np.exp(
            -0.5 * width * y**2 + 1.0j * p * y
        )
        second = ((-width * y + 1.0j * p) ** 2 - width) * gaussian
        hpsi += (-0.5 * second / model.mass_au)[:, None] * coefficients[None, :]
    hpsi += (
        psi @ model.H0.T
        + x[:, None] * (psi @ model.H1.T)
        + x[:, None] ** 2 * (psi @ model.H2.T)
    )

    tangents = []
    for mu in range(len(parameters)):
        step = 1.0e-6 * max(1.0, abs(parameters[mu]))
        plus = parameters.copy()
        minus = parameters.copy()
        plus[mu] += step
        minus[mu] -= step
        tangents.append((wavefunction(plus) - wavefunction(minus)) / (2.0 * step))
    numerical_metric = np.zeros_like(analytic.metric)
    numerical_rhs = np.zeros_like(analytic.rhs)
    for mu, left in enumerate(tangents):
        numerical_rhs[mu] = np.imag(
            np.trapezoid(np.sum(np.conj(left) * hpsi, axis=1), x)
        )
        for nu, right in enumerate(tangents):
            numerical_metric[mu, nu] = np.real(
                np.trapezoid(np.sum(np.conj(left) * right, axis=1), x)
            )
    assert np.allclose(analytic.metric, numerical_metric, atol=2.0e-9)
    assert np.allclose(analytic.rhs, numerical_rhs, atol=2.0e-9)


def test_single_packet_harmonic_tdvp_velocity_reduces_to_canonical_equations():
    mass = 900.0
    omega = 0.003
    force_constant = mass * omega**2
    model = QuadraticSpinHamiltonianV251(
        mass, np.zeros((1, 1)), np.zeros((1, 1)), [[0.5 * force_constant]]
    ).validate()
    state = FrozenGaussianSpinorStateV251(
        q=[-0.4],
        p=[2.3],
        widths=[mass * omega],
        coefficients=[[1.0 + 0.0j]],
    ).normalized()
    system = build_variational_metric_system_v251(state, model)
    assert system.solve_receipt.rank == state.parameter_count
    assert system.velocity[-2] == pytest.approx(state.p[0] / mass, abs=3.0e-15)
    assert system.velocity[-1] == pytest.approx(
        -force_constant * state.q[0], abs=3.0e-15
    )


def test_two_packet_metric_is_coupled_full_rank_and_residual_certified():
    state = _state()
    system = build_variational_metric_system_v251(state, _model())
    assert system.metric.shape == (state.parameter_count, state.parameter_count)
    assert system.solve_receipt.rank == state.parameter_count
    assert system.solve_receipt.nullity == 0
    assert system.solve_receipt.retained_condition_number < 600.0
    assert system.solve_receipt.linear_residual_relative < 2.0e-14
    coefficient_block = 2 * state.ngaussian * state.nstate
    assert np.max(np.abs(system.metric[:coefficient_block, coefficient_block:])) > 1.0e-3


def test_identical_packets_expose_compatible_metric_null_space():
    mass = 900.0
    omega = 0.003
    model = QuadraticSpinHamiltonianV251(
        mass,
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        [[0.5 * mass * omega**2]],
    ).validate()
    state = FrozenGaussianSpinorStateV251(
        q=[-0.4, -0.4],
        p=[2.3, 2.3],
        widths=[mass * omega, mass * omega],
        coefficients=[[0.5 + 0.1j], [0.5 - 0.1j]],
    ).normalized()
    system = build_variational_metric_system_v251(state, model)
    assert system.solve_receipt.rank == 4
    assert system.solve_receipt.nullity == 4
    assert system.solve_receipt.null_rhs_relative < 2.0e-15
    assert system.solve_receipt.linear_residual_relative < 2.0e-15


def test_incompatible_null_rhs_and_indefinite_metric_fail_closed():
    with pytest.raises(ValueError, match="incompatible with its null space"):
        solve_variational_metric_v251(
            np.diag([1.0, 0.0]), np.asarray([0.0, 1.0])
        )
    with pytest.raises(ValueError, match="positive semidefinite"):
        solve_variational_metric_v251(
            np.diag([1.0, -0.1]), np.asarray([0.0, 0.0])
        )


def test_implicit_midpoint_receipt_binds_endpoint_metric_and_solver():
    receipt = implicit_midpoint_tdvp_step_v251(_state(), _model(), 0.1)
    assert receipt.nonlinear_success is True
    assert receipt.nonlinear_status > 0
    assert receipt.nonlinear_function_evaluations >= 20
    assert receipt.nonlinear_residual_norm < 1.0e-14
    assert abs(receipt.norm_change) < 1.0e-13
    assert abs(receipt.energy_change_hartree) < 1.0e-12
    assert receipt.midpoint_system.solve_receipt.rank == 20
    assert receipt.as_dict()["model_fingerprint"] == _model().fingerprint()


def test_signed_implicit_midpoint_reversal_returns_the_initial_state():
    initial = _state()
    trajectory = run_frozen_width_multigaussian_tdvp_v251(
        initial, _model(), dt_au=0.1, steps=5
    )
    reverse = reverse_frozen_width_multigaussian_tdvp_v251(trajectory)
    assert _maximum_state_error(reverse.final_state, initial) < 2.0e-13
    assert trajectory.maximum_norm_drift < 5.0e-14
    assert trajectory.maximum_absolute_energy_drift_hartree < 5.0e-13
    assert trajectory.maximum_nonlinear_residual < 2.0e-14


def test_gaussian_permutation_covariance_is_exact_to_solver_precision():
    initial = _state()
    model = _model()
    base = run_frozen_width_multigaussian_tdvp_v251(
        initial, model, dt_au=0.1, steps=3
    )
    order = np.asarray([1, 0])
    permuted = run_frozen_width_multigaussian_tdvp_v251(
        initial.permuted(order), model, dt_au=0.1, steps=3
    )
    expected = base.final_state.permuted(order)
    assert _maximum_state_error(permuted.final_state, expected) < 2.0e-13


def test_constant_complex_electronic_gauge_covariance_is_complete():
    initial = _state()
    model = _model()
    unitary = random_unitary_v21(4, 25101)
    base = run_frozen_width_multigaussian_tdvp_v251(
        initial, model, dt_au=0.1, steps=3
    )
    transformed = run_frozen_width_multigaussian_tdvp_v251(
        initial.gauge_transformed(unitary),
        model.gauge_transformed(unitary),
        dt_au=0.1,
        steps=3,
    )
    expected = base.final_state.gauge_transformed(unitary)
    assert _maximum_state_error(transformed.final_state, expected) < 3.0e-13


def test_zero_soc_enabled_and_disabled_quadratic_paths_are_identical():
    enabled = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
        )
    )
    disabled = quadratic_spin_hamiltonian_from_provider_v251(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
        )
    )
    assert np.array_equal(enabled.H0, disabled.H0)
    assert np.array_equal(enabled.H1, disabled.H1)
    assert np.array_equal(enabled.H2, disabled.H2)
    left = run_frozen_width_multigaussian_tdvp_v251(
        _state(), enabled, dt_au=0.1, steps=2
    )
    right = run_frozen_width_multigaussian_tdvp_v251(
        _state(), disabled, dt_au=0.1, steps=2
    )
    assert _maximum_state_error(left.final_state, right.final_state) == 0.0


def test_timestep_refinement_has_the_implicit_midpoint_second_order_plateau():
    initial = _state()
    model = _model()
    trajectories = [
        run_frozen_width_multigaussian_tdvp_v251(
            initial, model, dt_au=dt, steps=int(round(0.4 / dt))
        )
        for dt in (0.1, 0.05, 0.025, 0.0125)
    ]
    changes = [
        _maximum_state_error(left.final_state, right.final_state)
        for left, right in zip(trajectories[:-1], trajectories[1:])
    ]
    ratios = [right / left for left, right in zip(changes[:-1], changes[1:])]
    assert all(0.245 < ratio < 0.255 for ratio in ratios)


def test_step_receipt_tampering_and_width_change_are_rejected():
    receipt = implicit_midpoint_tdvp_step_v251(_state(), _model(), 0.1)
    with pytest.raises(ValueError, match="nonlinear midpoint residual"):
        replace(
            receipt,
            nonlinear_residual=receipt.nonlinear_residual + 1.0e-3,
        ).validate()
    with pytest.raises(ValueError, match="function-evaluation count"):
        replace(receipt, nonlinear_function_evaluations=3.5).validate()
    broken_metric = replace(
        receipt.midpoint_system,
        metric=receipt.midpoint_system.metric + 1.0e-3 * np.eye(20),
    )
    with pytest.raises(ValueError):
        replace(receipt, midpoint_system=broken_metric).validate()
    with pytest.raises(ValueError, match="widths changed"):
        replace(
            receipt,
            end=replace(receipt.end, widths=receipt.end.widths + 1.0e-3),
        ).validate()


def test_nonlinear_nonconvergence_fails_closed():
    settings = replace(
        VariationalMetricSettingsV251(),
        nonlinear_max_function_evaluations=1,
    ).validate()
    with pytest.raises(RuntimeError, match="implicit midpoint TDVP solve failed"):
        implicit_midpoint_tdvp_step_v251(_state(), _model(), 0.8, settings=settings)
