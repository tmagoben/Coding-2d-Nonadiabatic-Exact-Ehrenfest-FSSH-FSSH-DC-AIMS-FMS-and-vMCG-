from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gaussian_dynamics.adaptive_multigaussian_tdvp_v252 import (
    ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
    VARIATIONAL_INTEGRATOR_V252,
    VARIATIONAL_METRIC_SOLVER_V252,
    WIDTH_COORDINATES_V252,
    AdaptiveVariationalSettingsV252,
    QuadraticSpinHamiltonianV251,
    ThawedGaussianSpinorStateV252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    adaptive_variational_energy_v252,
    build_adaptive_gaussian_spinor_matrices_v252,
    build_adaptive_variational_metric_system_v252,
    pack_adaptive_variational_parameters_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
    reverse_adaptive_width_multigaussian_tdvp_v252,
    run_adaptive_width_multigaussian_tdvp_v252,
    state_from_adaptive_variational_parameters_v252,
)
from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
    DoubletSOCConfigV220,
)
from gaussian_dynamics.complex_gauge_v21 import random_unitary_v21
from gaussian_dynamics.multigaussian_tdvp_v251 import (
    FrozenGaussianSpinorStateV251,
    run_frozen_width_multigaussian_tdvp_v251,
    solve_variational_metric_v251,
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
    return ThawedGaussianSpinorStateV252(
        q=[-0.65, 0.75],
        p=[5.0, -3.0],
        widths=[2.6, 2.1],
        chirps=[0.12, -0.08],
        coefficients=_coefficients(),
    ).normalized()


def _model():
    return quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220()
    )


def _maximum_state_error(left, right):
    return max(
        float(np.max(np.abs(left.q - right.q))),
        float(np.max(np.abs(left.p - right.p))),
        float(np.max(np.abs(left.widths - right.widths))),
        float(np.max(np.abs(left.chirps - right.chirps))),
        float(np.max(np.abs(left.coefficients - right.coefficients))),
        abs(float(left.time_au) - float(right.time_au)),
    )


def _harmonic_model():
    mass = 900.0
    omega = 0.003
    force_constant = mass * omega**2
    return QuadraticSpinHamiltonianV251(
        mass,
        np.zeros((1, 1)),
        np.zeros((1, 1)),
        [[0.5 * force_constant]],
    ).validate()


def test_scope_freezes_log_width_chirp_svd_and_implicit_midpoint():
    settings = AdaptiveVariationalSettingsV252().validate()
    assert "logarithmic widths" in ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252
    assert "chirp" in WIDTH_COORDINATES_V252
    assert settings.adaptive_gaussian_widths is True
    assert settings.integrator == VARIATIONAL_INTEGRATOR_V252
    assert settings.metric_solver == VARIATIONAL_METRIC_SOLVER_V252

    for field_name in (
        "spawning",
        "pruning",
        "coordinate_dependent_electronic_frame",
        "multidimensional_nuclear_motion",
        "full_width_matrices",
        "real_molecular_soc_provider",
    ):
        with pytest.raises(ValueError, match="does not admit"):
            replace(settings, **{field_name: True}).validate()
    with pytest.raises(ValueError, match="adaptive Gaussian widths as enabled"):
        replace(settings, adaptive_gaussian_widths=False).validate()
    with pytest.raises(ValueError, match="SVD metric solver is frozen"):
        replace(settings, metric_solver="normal equations").validate()
    with pytest.raises(ValueError, match="implicit integrator is frozen"):
        replace(settings, integrator="velocity Verlet").validate()


@pytest.mark.parametrize(
    "provider",
    [AnalyticSingletTripletSOCProviderV220(), AnalyticDoubletSOCProviderV220()],
)
def test_quadratic_intake_retains_complete_even_and_odd_soc(provider):
    model = quadratic_spin_hamiltonian_from_provider_v252(provider)
    assert model.nstate == 4
    assert model.physical_soc is True
    assert model.complete_spin_manifold is True
    assert model.source["fixed_frame_verified"] is True


def test_static_provider_intake_remains_closed():
    static = SimpleNamespace(evaluate_snapshot=lambda q: SimpleNamespace(matrices=1))
    with pytest.raises(TypeError, match="explicit operator provenance"):
        quadratic_spin_hamiltonian_from_provider_v252(static)


def test_state_log_width_pack_roundtrip_and_positivity_are_structural():
    state = _state()
    packed = pack_adaptive_variational_parameters_v252(state)
    restored = state_from_adaptive_variational_parameters_v252(
        packed,
        ngaussian=state.ngaussian,
        nstate=state.nstate,
        time_au=state.time_au,
    )
    assert _maximum_state_error(restored, state) < 5.0e-16
    width_start = 2 * state.ngaussian * state.nstate + 2 * state.ngaussian
    assert np.array_equal(
        packed[width_start : width_start + state.ngaussian], state.log_widths
    )
    with pytest.raises(ValueError, match="widths must be positive"):
        replace(state, widths=np.asarray([2.6, -2.1])).validate()
    with pytest.raises(ValueError, match="floating-point range"):
        state_from_adaptive_variational_parameters_v252(
            packed + np.eye(1, len(packed), width_start).reshape(-1) * 800.0,
            ngaussian=state.ngaussian,
            nstate=state.nstate,
            time_au=0.0,
        )


def test_exact_chirped_matrices_match_independent_grid_quadrature():
    model = QuadraticSpinHamiltonianV251(
        mass_au=7.0,
        H0=np.asarray([[0.04, 0.01 + 0.02j], [0.01 - 0.02j, 0.07]]),
        H1=np.asarray([[0.02, -0.008j], [0.008j, -0.01]]),
        H2=np.asarray([[0.006, 0.002], [0.002, 0.009]]),
    ).validate()
    state = ThawedGaussianSpinorStateV252(
        q=[-0.45, 0.62],
        p=[0.38, -0.27],
        widths=[1.35, 1.75],
        chirps=[0.31, -0.22],
        coefficients=[[0.8, 0.1j], [0.2 - 0.1j, 0.3]],
    ).normalized()
    analytic_overlap, analytic_hamiltonian = (
        build_adaptive_gaussian_spinor_matrices_v252(state, model)
    )

    x = np.linspace(-9.0, 9.0, 60001)
    gaussians = []
    kinetic_gaussians = []
    for q, p, width, chirp in zip(
        state.q, state.p, state.widths, state.chirps
    ):
        y = x - q
        gaussian = (width / np.pi) ** 0.25 * np.exp(
            -0.5 * width * y**2
            + 0.5j * chirp * y**2
            + 1.0j * p * y
        )
        z = width - 1.0j * chirp
        second = ((-z * y + 1.0j * p) ** 2 - z) * gaussian
        gaussians.append(gaussian)
        kinetic_gaussians.append(-0.5 * second / model.mass_au)
    numerical_overlap = np.zeros_like(analytic_overlap)
    numerical_hamiltonian = np.zeros_like(analytic_hamiltonian)
    for i in range(state.ngaussian):
        si = slice(i * state.nstate, (i + 1) * state.nstate)
        for j in range(state.ngaussian):
            sj = slice(j * state.nstate, (j + 1) * state.nstate)
            weights = np.conj(gaussians[i]) * gaussians[j]
            nuclear_overlap = np.trapezoid(weights, x)
            numerical_overlap[si, sj] = nuclear_overlap * np.eye(state.nstate)
            kinetic = np.trapezoid(
                np.conj(gaussians[i]) * kinetic_gaussians[j], x
            )
            potential = np.zeros((state.nstate, state.nstate), dtype=complex)
            for a in range(state.nstate):
                for b in range(state.nstate):
                    potential[a, b] = np.trapezoid(
                        weights
                        * (
                            model.H0[a, b]
                            + x * model.H1[a, b]
                            + x**2 * model.H2[a, b]
                        ),
                        x,
                    )
            numerical_hamiltonian[si, sj] = (
                kinetic * np.eye(state.nstate) + potential
            )
    assert np.allclose(analytic_overlap, numerical_overlap, atol=2.0e-11)
    assert np.allclose(analytic_hamiltonian, numerical_hamiltonian, atol=2.0e-11)


def test_adaptive_metric_and_rhs_match_independent_grid_tangents():
    model = QuadraticSpinHamiltonianV251(
        mass_au=8.0,
        H0=np.asarray([[0.03, 0.012 + 0.008j], [0.012 - 0.008j, 0.06]]),
        H1=np.asarray([[0.015, -0.004j], [0.004j, -0.009]]),
        H2=np.asarray([[0.005, 0.0015], [0.0015, 0.008]]),
    ).validate()
    state = ThawedGaussianSpinorStateV252(
        q=[-0.42, 0.58],
        p=[0.31, -0.22],
        widths=[1.4, 1.8],
        chirps=[0.27, -0.18],
        coefficients=[[0.72 + 0.05j, 0.08 - 0.12j], [0.19j, 0.31 - 0.04j]],
    ).normalized()
    analytic = build_adaptive_variational_metric_system_v252(state, model)
    parameters = pack_adaptive_variational_parameters_v252(state)
    x = np.linspace(-9.0, 9.0, 40001)

    def wavefunction(theta):
        local = state_from_adaptive_variational_parameters_v252(
            theta,
            ngaussian=state.ngaussian,
            nstate=state.nstate,
            time_au=0.0,
        )
        value = np.zeros((len(x), state.nstate), dtype=complex)
        for q, p, width, chirp, coefficients in zip(
            local.q,
            local.p,
            local.widths,
            local.chirps,
            local.coefficients,
        ):
            y = x - q
            gaussian = (width / np.pi) ** 0.25 * np.exp(
                -0.5 * width * y**2
                + 0.5j * chirp * y**2
                + 1.0j * p * y
            )
            value += gaussian[:, None] * coefficients[None, :]
        return value

    psi = wavefunction(parameters)
    hpsi = np.zeros_like(psi)
    for q, p, width, chirp, coefficients in zip(
        state.q,
        state.p,
        state.widths,
        state.chirps,
        state.coefficients,
    ):
        y = x - q
        gaussian = (width / np.pi) ** 0.25 * np.exp(
            -0.5 * width * y**2
            + 0.5j * chirp * y**2
            + 1.0j * p * y
        )
        z = width - 1.0j * chirp
        second = ((-z * y + 1.0j * p) ** 2 - z) * gaussian
        hpsi += (-0.5 * second / model.mass_au)[:, None] * coefficients[None, :]
    hpsi += (
        psi @ model.H0.T
        + x[:, None] * (psi @ model.H1.T)
        + x[:, None] ** 2 * (psi @ model.H2.T)
    )

    tangents = []
    for mu in range(len(parameters)):
        step = 8.0e-7 * max(1.0, abs(parameters[mu]))
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
    assert np.allclose(analytic.metric, numerical_metric, atol=3.0e-9)
    assert np.allclose(analytic.rhs, numerical_rhs, atol=3.0e-9)


def test_single_packet_harmonic_velocity_matches_thawed_gaussian_equations():
    model = _harmonic_model()
    state = ThawedGaussianSpinorStateV252(
        q=[-0.4],
        p=[2.3],
        widths=[2.1],
        chirps=[0.37],
        coefficients=[[1.0]],
    ).normalized()
    system = build_adaptive_variational_metric_system_v252(state, model)
    mass = model.mass_au
    force_constant = 2.0 * model.H2[0, 0].real
    assert system.solve_receipt.rank == state.parameter_count
    assert system.velocity[-4] == pytest.approx(state.p[0] / mass, abs=3.0e-15)
    assert system.velocity[-3] == pytest.approx(
        -force_constant * state.q[0], abs=3.0e-15
    )
    assert system.velocity[-2] == pytest.approx(
        -2.0 * state.chirps[0] / mass, abs=3.0e-15
    )
    assert system.velocity[-1] == pytest.approx(
        (state.widths[0] ** 2 - state.chirps[0] ** 2) / mass
        - force_constant,
        abs=3.0e-15,
    )


def test_coherent_harmonic_state_reduces_to_v0251_frozen_trajectory():
    model = _harmonic_model()
    width = model.mass_au * 0.003
    adaptive = ThawedGaussianSpinorStateV252(
        [-0.4], [2.3], [width], [0.0], [[1.0]]
    ).normalized()
    frozen = FrozenGaussianSpinorStateV251(
        [-0.4], [2.3], [width], [[1.0]]
    ).normalized()
    adaptive_trajectory = run_adaptive_width_multigaussian_tdvp_v252(
        adaptive, model, dt_au=0.1, steps=4
    )
    frozen_trajectory = run_frozen_width_multigaussian_tdvp_v251(
        frozen, model, dt_au=0.1, steps=4
    )
    assert np.max(
        np.abs(adaptive_trajectory.final_state.q - frozen_trajectory.final_state.q)
    ) < 2.0e-14
    assert np.max(
        np.abs(adaptive_trajectory.final_state.p - frozen_trajectory.final_state.p)
    ) < 2.0e-14
    assert np.max(
        np.abs(
            adaptive_trajectory.final_state.coefficients
            - frozen_trajectory.final_state.coefficients
        )
    ) < 3.0e-14
    assert adaptive_trajectory.final_state.widths[0] == pytest.approx(width, abs=2.0e-14)
    assert adaptive_trajectory.final_state.chirps[0] == pytest.approx(0.0, abs=2.0e-14)


def test_two_packet_adaptive_metric_is_coupled_full_rank():
    state = _state()
    system = build_adaptive_variational_metric_system_v252(state, _model())
    assert system.metric.shape == (state.parameter_count, state.parameter_count)
    assert system.solve_receipt.rank == state.parameter_count
    assert system.solve_receipt.nullity == 0
    assert system.solve_receipt.retained_condition_number < 5000.0
    assert system.solve_receipt.linear_residual_relative < 2.0e-14
    coefficient_block = 2 * state.ngaussian * state.nstate
    assert np.max(np.abs(system.metric[:coefficient_block, coefficient_block:])) > 1.0e-3


def test_duplicate_adaptive_packets_expose_compatible_null_space():
    state = ThawedGaussianSpinorStateV252(
        q=[-0.2, -0.2],
        p=[0.7, 0.7],
        widths=[2.4, 2.4],
        chirps=[0.2, 0.2],
        coefficients=[[0.6 + 0.1j], [0.3 - 0.2j]],
    ).normalized()
    system = build_adaptive_variational_metric_system_v252(state, _harmonic_model())
    assert system.solve_receipt.rank == 6
    assert system.solve_receipt.nullity == 6
    assert system.solve_receipt.null_rhs_relative < 2.0e-15
    assert system.solve_receipt.linear_residual_relative < 2.0e-14


def test_incompatible_null_rhs_and_indefinite_metric_still_fail_closed():
    with pytest.raises(ValueError, match="incompatible with its null space"):
        solve_variational_metric_v251(
            np.diag([1.0, 0.0]),
            np.asarray([0.0, 1.0]),
            settings=AdaptiveVariationalSettingsV252(),
        )
    with pytest.raises(ValueError, match="positive semidefinite"):
        solve_variational_metric_v251(
            np.diag([1.0, -0.1]),
            np.zeros(2),
            settings=AdaptiveVariationalSettingsV252(),
        )


def test_adaptive_midpoint_receipt_binds_width_chirp_metric_and_endpoint():
    receipt = adaptive_implicit_midpoint_tdvp_step_v252(_state(), _model(), 0.05)
    assert receipt.nonlinear_success is True
    assert receipt.nonlinear_function_evaluations >= 24
    assert receipt.nonlinear_residual_norm < 2.0e-14
    assert abs(receipt.norm_change) < 2.0e-13
    assert abs(receipt.energy_change_hartree) < 2.0e-12
    assert receipt.midpoint_system.solve_receipt.rank == 24
    assert receipt.maximum_log_width_change < 2.0e-5
    assert receipt.as_dict()["model_fingerprint"] == _model().fingerprint()


def test_signed_adaptive_midpoint_reversal_returns_initial_state():
    initial = _state()
    trajectory = run_adaptive_width_multigaussian_tdvp_v252(
        initial, _model(), dt_au=0.05, steps=5
    )
    reverse = reverse_adaptive_width_multigaussian_tdvp_v252(trajectory)
    assert _maximum_state_error(reverse.final_state, initial) < 3.0e-13
    assert trajectory.maximum_norm_drift < 2.0e-13
    assert trajectory.maximum_absolute_energy_drift_hartree < 2.0e-12
    assert trajectory.maximum_nonlinear_residual < 2.0e-14


def test_gaussian_permutation_covariance_includes_width_and_chirp():
    initial = _state()
    base = run_adaptive_width_multigaussian_tdvp_v252(
        initial, _model(), dt_au=0.05, steps=3
    )
    order = np.asarray([1, 0])
    permuted = run_adaptive_width_multigaussian_tdvp_v252(
        initial.permuted(order), _model(), dt_au=0.05, steps=3
    )
    assert _maximum_state_error(
        permuted.final_state, base.final_state.permuted(order)
    ) < 3.0e-13


def test_constant_complex_electronic_gauge_covariance_is_complete():
    initial = _state()
    model = _model()
    unitary = random_unitary_v21(4, 25201)
    base = run_adaptive_width_multigaussian_tdvp_v252(
        initial, model, dt_au=0.05, steps=3
    )
    transformed = run_adaptive_width_multigaussian_tdvp_v252(
        initial.gauge_transformed(unitary),
        model.gauge_transformed(unitary),
        dt_au=0.05,
        steps=3,
    )
    assert _maximum_state_error(
        transformed.final_state, base.final_state.gauge_transformed(unitary)
    ) < 4.0e-13


def test_zero_soc_enabled_and_disabled_paths_are_identical():
    enabled = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=True)
        )
    )
    disabled = quadratic_spin_hamiltonian_from_provider_v252(
        AnalyticDoubletSOCProviderV220(
            DoubletSOCConfigV220(soc_scale=0.0, soc_enabled=False)
        )
    )
    left = run_adaptive_width_multigaussian_tdvp_v252(
        _state(), enabled, dt_au=0.05, steps=2
    )
    right = run_adaptive_width_multigaussian_tdvp_v252(
        _state(), disabled, dt_au=0.05, steps=2
    )
    assert _maximum_state_error(left.final_state, right.final_state) == 0.0


def test_adaptive_timestep_refinement_has_second_order_plateau():
    trajectories = [
        run_adaptive_width_multigaussian_tdvp_v252(
            _state(), _model(), dt_au=dt, steps=int(round(0.2 / dt))
        )
        for dt in (0.1, 0.05, 0.025, 0.0125)
    ]
    changes = [
        _maximum_state_error(left.final_state, right.final_state)
        for left, right in zip(trajectories[:-1], trajectories[1:])
    ]
    ratios = [right / left for left, right in zip(changes[:-1], changes[1:])]
    assert all(0.245 < ratio < 0.255 for ratio in ratios)


def test_receipt_tampering_and_width_domain_violations_are_rejected():
    receipt = adaptive_implicit_midpoint_tdvp_step_v252(_state(), _model(), 0.05)
    with pytest.raises(ValueError, match="adaptive nonlinear residual"):
        replace(
            receipt,
            nonlinear_residual=receipt.nonlinear_residual + 1.0e-3,
        ).validate()
    broken_metric = replace(
        receipt.midpoint_system,
        metric=receipt.midpoint_system.metric + 1.0e-3 * np.eye(24),
    )
    with pytest.raises(ValueError):
        replace(receipt, midpoint_system=broken_metric).validate()
    with pytest.raises(ValueError, match="configured minimum"):
        replace(
            receipt,
            end=replace(receipt.end, widths=np.asarray([1.0e-10, 2.1])),
        ).validate()
    assert np.isfinite(adaptive_variational_energy_v252(receipt.end, receipt.model))


def test_nonlinear_nonconvergence_fails_closed():
    settings = replace(
        AdaptiveVariationalSettingsV252(),
        nonlinear_max_function_evaluations=1,
    ).validate()
    with pytest.raises(
        RuntimeError, match="adaptive implicit midpoint TDVP solve failed"
    ):
        adaptive_implicit_midpoint_tdvp_step_v252(
            _state(), _model(), 0.8, settings=settings
        )
