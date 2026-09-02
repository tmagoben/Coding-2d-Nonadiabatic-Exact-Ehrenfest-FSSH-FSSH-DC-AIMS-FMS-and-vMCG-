"""Deterministic scientific validation for the v0.26.0 multidimensional release."""

from dataclasses import dataclass
import json

import numpy as np

from .adaptive_multigaussian_tdvp_v252 import (
    ThawedGaussianSpinorStateV252,
    build_adaptive_gaussian_spinor_matrices_v252,
    build_adaptive_variational_metric_system_v252,
)
from .multigaussian_tdvp_v251 import QuadraticSpinHamiltonianV251
from .multidimensional_basis_adaptation_v260 import (
    ControlledMultidimensionalBasisSettingsV260,
    MultidimensionalSpawnCandidateV260,
    V260_MULTIDIMENSIONAL_BASIS_CLAIMS,
    adapt_multidimensional_basis_once_v260,
    evaluate_multidimensional_spawn_candidate_v260,
    generate_multidimensional_spawn_candidates_v260,
    metric_compatible_activation_mask_v260,
    run_controlled_multidimensional_dynamics_v260,
)
from .multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    V260_MULTIDIMENSIONAL_TDVP_CLAIMS,
    build_multidimensional_gaussian_matrices_v260,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
    multidimensional_reduced_density_v260,
    multidimensional_state_on_grid_v260,
    multidimensional_variational_energy_v260,
    pack_multidimensional_parameters_v260,
    run_multidimensional_tdvp_v260,
    state_from_multidimensional_parameters_v260,
)
from .multidimensional_soc_v260 import (
    EXACT_GRID_SCHEMA_V260,
    V260_EXACT_GRID_CLAIMS,
    ExactGridSettingsV260,
    QuadraticSpinHamiltonianNDV260,
    UniformGrid2DV260,
    exact_grid_boundary_probability_v260,
    exact_grid_split_step_v260,
    initial_gaussian_spinor_2d_v260,
    kramers_doublet_ci_soc_model_v260,
    normalize_spinor_grid_v260,
    phase_aligned_grid_error_v260,
    run_exact_grid_ci_soc_v260,
    singlet_triplet_ci_soc_model_v260,
    two_state_ci_soc_model_v260,
)


MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260 = "gnd-multidimensional-ci-soc-validation-v0.26.0"


def _maximum_parameter_error_v260(left, right):
    return float(
        np.max(
            np.abs(
                pack_multidimensional_parameters_v260(left)
                - pack_multidimensional_parameters_v260(right)
            )
        )
    )


def _random_unitary_v260(dimension, seed):
    rng = np.random.default_rng(int(seed))
    matrix = rng.normal(size=(dimension, dimension)) + 1.0j * rng.normal(
        size=(dimension, dimension)
    )
    q, r = np.linalg.qr(matrix)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0.0, phases / np.abs(phases), 1.0)
    return q @ np.diag(phases.conj())


def _velocity_transform_error_v260(state, model, transform_state, transform_model):
    system = build_multidimensional_metric_system_v260(state, model)
    transformed_state = transform_state(state)
    transformed_system = build_multidimensional_metric_system_v260(
        transformed_state, transform_model(model)
    )
    epsilon = 1.0e-6
    displaced = state_from_multidimensional_parameters_v260(
        pack_multidimensional_parameters_v260(state) + epsilon * system.velocity,
        ngaussian=state.ngaussian,
        ndim=state.ndim,
        nstate=state.nstate,
        time_au=state.time_au,
    )
    expected = (
        pack_multidimensional_parameters_v260(transform_state(displaced))
        - pack_multidimensional_parameters_v260(transformed_state)
    ) / epsilon
    return float(np.max(np.abs(expected - transformed_system.velocity)))


def _independent_grid_matrix_error_v260(state, model):
    """Compare analytic S/H against an FFT-grid quadrature with no moment code."""

    grid = UniformGrid2DV260.from_bounds((-9.0, 9.0), (-9.0, 9.0), (160, 160))
    points = grid.mesh()
    gaussians = []
    for packet in range(state.ngaussian):
        displacement = points - state.q[packet]
        normalization = float(np.prod(state.widths[packet] / np.pi) ** 0.25)
        exponent = np.sum(
            -0.5 * state.widths[packet] * displacement**2
            + 0.5j * state.chirps[packet] * displacement**2
            + 1.0j * state.p[packet] * displacement,
            axis=-1,
        )
        gaussians.append(normalization * np.exp(exponent))
    gaussians = np.asarray(gaussians)
    potential = model.hamiltonian(points)
    kx = 2.0 * np.pi * np.fft.fftfreq(len(grid.x), d=grid.dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(len(grid.y), d=grid.dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    wavevectors = np.stack((KX, KY), axis=-1)
    kinetic_energy = 0.5 * np.einsum(
        "...a,ab,...b->...",
        wavevectors,
        model.inverse_mass_matrix_au,
        wavevectors,
        optimize=True,
    )
    dimension = state.ngaussian * state.nstate
    overlap_grid = np.zeros((dimension, dimension), dtype=complex)
    hamiltonian_grid = np.zeros((dimension, dimension), dtype=complex)
    for packet_j in range(state.ngaussian):
        kinetic = np.fft.ifftn(
            kinetic_energy * np.fft.fftn(gaussians[packet_j])
        )
        for electronic_j in range(state.nstate):
            ket_index = packet_j * state.nstate + electronic_j
            action = potential[..., :, electronic_j] * gaussians[packet_j][..., None]
            action[..., electronic_j] += kinetic
            for packet_i in range(state.ngaussian):
                for electronic_i in range(state.nstate):
                    bra_index = packet_i * state.nstate + electronic_i
                    overlap_grid[bra_index, ket_index] = (
                        np.vdot(gaussians[packet_i], gaussians[packet_j])
                        * grid.volume_element
                        if electronic_i == electronic_j
                        else 0.0
                    )
                    hamiltonian_grid[bra_index, ket_index] = (
                        np.vdot(gaussians[packet_i], action[..., electronic_i])
                        * grid.volume_element
                    )
    overlap_analytic, hamiltonian_analytic = build_multidimensional_gaussian_matrices_v260(
        state, model
    )
    return (
        float(np.max(np.abs(overlap_grid - overlap_analytic))),
        float(np.max(np.abs(hamiltonian_grid - hamiltonian_analytic))),
    )


@dataclass(frozen=True)
class MultidimensionalValidationEvidenceV260:
    metrics: dict
    thresholds: dict
    checks: dict
    exact_grid_fingerprint: str
    controlled_trajectory_fingerprint: str

    @property
    def passed(self):
        return all(self.checks.values())

    @property
    def check_count(self):
        return len(self.checks)

    def validate(self):
        if not self.checks:
            raise ValueError("v0.26.0 validation evidence requires checks.")
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.26.0 validation gate must be a native Boolean.")
        if not self.passed:
            failed = ", ".join(name for name, value in self.checks.items() if not value)
            raise ValueError("v0.26.0 validation failed: " + failed)
        for mapping in (self.metrics, self.thresholds):
            for name, value in mapping.items():
                if not np.isfinite(float(value)):
                    raise ValueError(f"validation scalar {name!r} is non-finite.")
        if len(self.exact_grid_fingerprint) != 64 or len(self.controlled_trajectory_fingerprint) != 64:
            raise ValueError("validation evidence fingerprints must be SHA-256 strings.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260,
            "metrics": {name: float(value) for name, value in self.metrics.items()},
            "thresholds": {name: float(value) for name, value in self.thresholds.items()},
            "checks": dict(self.checks),
            "passed": self.passed,
            "check_count": self.check_count,
            "fingerprints": {
                "exact_grid": self.exact_grid_fingerprint,
                "controlled_trajectory": self.controlled_trajectory_fingerprint,
            },
            "claims": {
                "exact_grid": dict(V260_EXACT_GRID_CLAIMS),
                "tdvp": dict(V260_MULTIDIMENSIONAL_TDVP_CLAIMS),
                "basis": dict(V260_MULTIDIMENSIONAL_BASIS_CLAIMS),
            },
        }

    def fingerprint(self):
        payload = json.dumps(
            self.as_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
        import hashlib

        return hashlib.sha256(payload).hexdigest()


def run_multidimensional_validation_evidence_v260():
    # --- Model algebra, complete spin spaces, and symmetry contracts.
    model = two_state_ci_soc_model_v260(
        mass_au=(50.0, 50.0),
        kappa=0.04,
        coupling=0.04,
        frequencies=(0.03, 0.03),
        soc_scale=0.01,
    )
    zero_soc_model = two_state_ci_soc_model_v260(
        mass_au=(50.0, 50.0),
        kappa=0.04,
        coupling=0.04,
        frequencies=(0.03, 0.03),
        soc_scale=0.0,
    )
    doublet_model = kramers_doublet_ci_soc_model_v260()
    singlet_triplet_model = singlet_triplet_ci_soc_model_v260()
    origin_gap = float(np.diff(np.linalg.eigvalsh(model.H0))[0])
    zero_origin_gap = float(np.diff(np.linalg.eigvalsh(zero_soc_model.H0))[0])
    coordinate = np.asarray([0.37, -0.21])
    epsilon = 1.0e-6
    finite_difference = np.asarray(
        [
            (
                model.hamiltonian(coordinate + epsilon * np.eye(2)[axis])
                - model.hamiltonian(coordinate - epsilon * np.eye(2)[axis])
            )
            / (2.0 * epsilon)
            for axis in range(2)
        ]
    )
    derivative_error = float(np.max(np.abs(finite_difference - model.derivative(coordinate))))
    angle = 0.371
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    rotated_model = model.coordinate_rotated(rotation)
    coordinate_new = np.asarray([-0.13, 0.44])
    coordinate_transform_error = float(
        np.max(
            np.abs(
                rotated_model.hamiltonian(coordinate_new)
                - model.hamiltonian(coordinate_new @ rotation)
            )
        )
    )
    time_reversal_unitary = np.kron(
        np.eye(2), np.asarray([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
    )
    doublet_hamiltonian = doublet_model.hamiltonian(coordinate)
    doublet_time_reversal_error = float(
        np.max(
            np.abs(
                time_reversal_unitary
                @ doublet_hamiltonian.conj()
                @ time_reversal_unitary.conj().T
                - doublet_hamiltonian
            )
        )
    )
    doublet_energies = np.linalg.eigvalsh(doublet_hamiltonian)
    kramers_splitting = float(
        max(abs(doublet_energies[1] - doublet_energies[0]), abs(doublet_energies[3] - doublet_energies[2]))
    )

    # --- Independent exact-grid oracle, convergence, reversibility, and boundaries.
    grid = UniformGrid2DV260.from_bounds((-6.0, 6.0), (-6.0, 6.0), (64, 64))
    psi0 = initial_gaussian_spinor_2d_v260(
        grid,
        [1.0, 0.0],
        center=(-0.25, 0.0),
        momentum=(3.0, 0.0),
        widths=(2.0, 2.0),
    )
    exact = run_exact_grid_ci_soc_v260(
        model,
        grid,
        psi0,
        settings=ExactGridSettingsV260(dt_au=0.01, steps=6, store_every=6),
    )
    reverse = exact.final_state.copy()
    for _ in range(6):
        reverse = exact_grid_split_step_v260(reverse, model, grid, -0.01)
    grid_reversal_error = phase_aligned_grid_error_v260(psi0, reverse, grid)
    boundary_probability = max(
        exact_grid_boundary_probability_v260(psi0, grid),
        exact_grid_boundary_probability_v260(exact.final_state, grid),
    )
    refinement_runs = []
    for dt, steps in ((0.04, 5), (0.02, 10), (0.01, 20), (0.005, 40)):
        refinement_runs.append(
            run_exact_grid_ci_soc_v260(
                model,
                grid,
                psi0,
                settings=ExactGridSettingsV260(dt_au=dt, steps=steps, store_every=steps),
            )
        )
    refinement_reference = refinement_runs[-1].final_state
    refinement_errors = [
        phase_aligned_grid_error_v260(refinement_reference, item.final_state, grid)
        for item in refinement_runs[:-1]
    ]
    grid_order_coarse = float(np.log2(refinement_errors[0] / refinement_errors[1]))
    grid_order_fine = float(np.log2(refinement_errors[1] / refinement_errors[2]))
    doublet_grid = UniformGrid2DV260.from_bounds((-5.0, 5.0), (-5.0, 5.0), (24, 24))
    doublet_exact = run_exact_grid_ci_soc_v260(
        doublet_model,
        doublet_grid,
        initial_gaussian_spinor_2d_v260(doublet_grid, [1.0, 0.0, 0.0, 0.0]),
        settings=ExactGridSettingsV260(dt_au=0.01, steps=2, store_every=2),
    )
    st_grid = UniformGrid2DV260.from_bounds((-5.0, 5.0), (-5.0, 5.0), (24, 24))
    singlet_triplet_exact = run_exact_grid_ci_soc_v260(
        singlet_triplet_model,
        st_grid,
        initial_gaussian_spinor_2d_v260(st_grid, [1.0, 0.0, 0.0, 0.0, 0.0]),
        settings=ExactGridSettingsV260(dt_au=0.01, steps=2, store_every=2),
    )

    # --- Analytic multidimensional Gaussian algebra and inherited 1D reduction.
    quadrature_state = DiagonalGaussianSpinorStateV260(
        q=[[-0.7, 0.2], [0.8, -0.3]],
        p=[[0.5, -0.2], [-0.3, 0.4]],
        widths=[[1.1, 0.9], [0.8, 1.2]],
        chirps=[[0.1, -0.04], [-0.05, 0.08]],
        coefficients=[[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    quadrature_model = two_state_ci_soc_model_v260(
        mass_au=(900.0, 700.0),
        kappa=0.003,
        coupling=0.004,
        frequencies=(0.02, 0.025),
        soc_scale=0.002,
    )
    overlap_quadrature_error, hamiltonian_quadrature_error = _independent_grid_matrix_error_v260(
        quadrature_state, quadrature_model
    )
    H0 = np.asarray([[0.01, 0.002j], [-0.002j, -0.01]])
    H1 = np.asarray([[0.003, 0.004], [0.004, -0.002]])
    H2 = np.asarray([[0.0005, 0.0], [0.0, 0.0007]])
    model_1d_old = QuadraticSpinHamiltonianV251(900.0, H0, H1, H2).validate()
    model_1d_new = QuadraticSpinHamiltonianNDV260(
        [[900.0]], H0, H1[None, :, :], H2[None, None, :, :]
    ).validate()
    state_1d_old = ThawedGaussianSpinorStateV252(
        q=[-0.7, 0.8],
        p=[0.5, -0.3],
        widths=[1.1, 0.8],
        chirps=[0.1, -0.05],
        coefficients=[[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    state_1d_new = DiagonalGaussianSpinorStateV260(
        state_1d_old.q[:, None],
        state_1d_old.p[:, None],
        state_1d_old.widths[:, None],
        state_1d_old.chirps[:, None],
        state_1d_old.coefficients,
    ).validate(require_normalized=True)
    old_matrices = build_adaptive_gaussian_spinor_matrices_v252(state_1d_old, model_1d_old)
    new_matrices = build_multidimensional_gaussian_matrices_v260(state_1d_new, model_1d_new)
    old_metric = build_adaptive_variational_metric_system_v252(state_1d_old, model_1d_old)
    new_metric = build_multidimensional_metric_system_v260(state_1d_new, model_1d_new)
    reduction_overlap_error = float(np.max(np.abs(old_matrices[0] - new_matrices[0])))
    reduction_hamiltonian_error = float(np.max(np.abs(old_matrices[1] - new_matrices[1])))
    reduction_metric_error = float(np.max(np.abs(old_metric.metric - new_metric.metric)))
    reduction_rhs_error = float(np.max(np.abs(old_metric.rhs - new_metric.rhs)))
    reduction_velocity_error = float(np.max(np.abs(old_metric.velocity - new_metric.velocity)))

    # --- TDVP symmetry, midpoint, and exact-reference comparison.
    initial = DiagonalGaussianSpinorStateV260(
        q=[[-0.25, 0.0]],
        p=[[3.0, 0.0]],
        widths=[[2.0, 2.0]],
        chirps=[[0.0, 0.0]],
        coefficients=[[1.0, 0.0]],
    ).normalized()
    metric_system = build_multidimensional_metric_system_v260(initial, model)
    metric_minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(metric_system.metric)))
    forward = multidimensional_implicit_midpoint_step_v260(initial, model, 0.01)
    backward = multidimensional_implicit_midpoint_step_v260(forward.end, model, -0.01)
    tdvp_reversal_error = _maximum_parameter_error_v260(initial, backward.end)
    unitary = _random_unitary_v260(2, 26001)
    gauge_velocity_error = _velocity_transform_error_v260(
        initial,
        model,
        lambda state: state.gauge_transformed(unitary),
        lambda item: item.gauge_transformed(unitary),
    )
    permutation_state = quadrature_state
    order = np.asarray([1, 0])
    permutation_velocity_error = _velocity_transform_error_v260(
        permutation_state,
        quadrature_model,
        lambda state: state.permuted(order),
        lambda item: item,
    )
    coordinate_swap = np.asarray([[0.0, 1.0], [-1.0, 0.0]])
    coordinate_velocity_error = _velocity_transform_error_v260(
        quadrature_state,
        quadrature_model,
        lambda state: state.coordinate_rotated(coordinate_swap),
        lambda item: item.coordinate_rotated(coordinate_swap),
    )
    one_packet = run_multidimensional_tdvp_v260(initial, model, 0.01, 6)
    controlled = run_controlled_multidimensional_dynamics_v260(initial, model, 0.01, 6)
    one_packet_grid = normalize_spinor_grid_v260(
        multidimensional_state_on_grid_v260(one_packet.final_state, grid), grid
    )
    controlled_grid = normalize_spinor_grid_v260(
        multidimensional_state_on_grid_v260(controlled.final_state, grid), grid
    )
    one_packet_exact_error = phase_aligned_grid_error_v260(
        exact.final_state, one_packet_grid, grid
    )
    controlled_exact_error = phase_aligned_grid_error_v260(
        exact.final_state, controlled_grid, grid
    )
    one_packet_density_error = float(
        np.linalg.norm(
            exact.reduced_densities[-1]
            - multidimensional_reduced_density_v260(one_packet.final_state)
        )
    )
    controlled_density_error = float(
        np.linalg.norm(
            exact.reduced_densities[-1]
            - multidimensional_reduced_density_v260(controlled.final_state)
        )
    )

    # --- Lifecycle gates: duplicate rejection, prune, merge, activation, reductions.
    candidates = generate_multidimensional_spawn_candidates_v260(initial)
    candidate_evaluations = [
        evaluate_multidimensional_spawn_candidate_v260(initial, model, candidate)
        for candidate in candidates
    ]
    duplicate_candidate = MultidimensionalSpawnCandidateV260(
        initial.q[0], initial.p[0], initial.widths[0], initial.chirps[0], 0, "position", 0, 1
    )
    duplicate_evaluation = evaluate_multidimensional_spawn_candidate_v260(
        initial, model, duplicate_candidate
    )
    prune_state = DiagonalGaussianSpinorStateV260(
        q=[[-2.0, 0.0], [2.0, 0.0]],
        p=[[0.0, 0.0], [0.0, 0.0]],
        widths=[[2.0, 2.0], [2.0, 2.0]],
        chirps=[[0.0, 0.0], [0.0, 0.0]],
        coefficients=[[1.0, 0.0], [1.0e-7, 0.0]],
    ).normalized()
    prune_event = adapt_multidimensional_basis_once_v260(
        prune_state,
        zero_soc_model,
        packet_ids=("g000000", "g000001"),
        packet_ages=(64, 64),
        next_packet_serial=2,
    )
    merge_state = DiagonalGaussianSpinorStateV260(
        q=[[0.0, 0.0], [1.0e-4, 0.0]],
        p=[[0.0, 0.0], [0.0, 0.0]],
        widths=[[2.0, 2.0], [2.0, 2.0]],
        chirps=[[0.0, 0.0], [0.0, 0.0]],
        coefficients=[[0.7, 0.0], [0.3, 0.0]],
    ).normalized()
    merge_event = adapt_multidimensional_basis_once_v260(
        merge_state,
        zero_soc_model,
        packet_ids=("g000000", "g000001"),
        packet_ages=(2, 2),
        next_packet_serial=2,
    )
    spawn_event = adapt_multidimensional_basis_once_v260(initial, model)
    dormant_state = spawn_event.after
    dormant_mask = metric_compatible_activation_mask_v260(
        dormant_state, model, locked_active_mask=[True, False]
    )
    activated_coefficients = dormant_state.coefficients.copy()
    activated_coefficients[-1, 1] = 0.01
    activation_state = DiagonalGaussianSpinorStateV260(
        dormant_state.q,
        dormant_state.p,
        dormant_state.widths,
        dormant_state.chirps,
        activated_coefficients,
    ).normalized()
    activated_mask = metric_compatible_activation_mask_v260(
        activation_state, model, locked_active_mask=[True, False]
    )
    inactive_step = multidimensional_implicit_midpoint_step_v260(
        dormant_state, model, 0.01, active_shape_mask=dormant_mask
    )
    dormant_shape_drift = max(
        float(np.max(np.abs(inactive_step.end.q[-1] - inactive_step.start.q[-1]))),
        float(np.max(np.abs(inactive_step.end.p[-1] - inactive_step.start.p[-1]))),
        float(np.max(np.abs(inactive_step.end.widths[-1] - inactive_step.start.widths[-1]))),
        float(np.max(np.abs(inactive_step.end.chirps[-1] - inactive_step.start.chirps[-1]))),
    )
    no_event_settings = ControlledMultidimensionalBasisSettingsV260(
        spawn_residual_capture_threshold=1.0e3,
        adapt_every_steps=1,
    )
    no_event = run_controlled_multidimensional_dynamics_v260(
        initial, model, 0.01, 1, settings=no_event_settings
    )
    direct = multidimensional_implicit_midpoint_step_v260(initial, model, 0.01)
    no_event_reduction_error = _maximum_parameter_error_v260(no_event.final_state, direct.end)
    zero_enabled = adapt_multidimensional_basis_once_v260(initial, zero_soc_model)
    zero_disabled = adapt_multidimensional_basis_once_v260(
        initial,
        two_state_ci_soc_model_v260(
            mass_au=(50.0, 50.0),
            kappa=0.04,
            coupling=0.04,
            frequencies=(0.03, 0.03),
            soc_scale=-0.0,
        ),
    )
    zero_soc_score_error = float(
        np.max(
            np.abs(
                np.asarray([item.residual_capture for item in zero_enabled.candidate_evaluations])
                - np.asarray([item.residual_capture for item in zero_disabled.candidate_evaluations])
            )
        )
    )

    thresholds = {
        "model_derivative_error": 2.0e-10,
        "coordinate_transform_error": 2.0e-11,
        "time_reversal_symmetry_error": 2.0e-11,
        "kramers_splitting": 2.0e-11,
        "grid_norm_drift": 2.0e-12,
        "grid_reversal_error": 2.0e-11,
        "grid_boundary_probability": 1.0e-12,
        "minimum_grid_order": 1.8,
        "grid_energy_drift": 2.0e-10,
        "quadrature_overlap_error": 2.0e-10,
        "quadrature_hamiltonian_error": 2.0e-9,
        "one_dimensional_reduction_error": 3.0e-13,
        "metric_minimum_eigenvalue": -3.0e-10,
        "metric_linear_residual": 3.0e-9,
        "tdvp_reversal_error": 2.0e-8,
        "covariance_error": 2.0e-8,
        "tdvp_norm_drift": 1.0e-8,
        "tdvp_energy_drift": 1.0e-8,
        "maximum_exact_wavefunction_error": 3.0e-3,
        "maximum_density_error": 2.0e-5,
        "adaptive_improvement_ratio": 0.8,
        "projection_loss": 2.0e-7,
        "projection_energy_jump": 2.0e-6,
        "dormant_shape_drift": 2.0e-13,
        "no_event_reduction_error": 3.0e-13,
        "zero_soc_error": 3.0e-13,
    }
    metrics = {
        "origin_soc_gap_hartree": origin_gap,
        "zero_soc_origin_gap_hartree": zero_origin_gap,
        "model_derivative_error": derivative_error,
        "coordinate_transform_error": coordinate_transform_error,
        "doublet_time_reversal_error": doublet_time_reversal_error,
        "maximum_kramers_splitting_hartree": kramers_splitting,
        "exact_grid_norm_drift": exact.maximum_norm_drift,
        "exact_grid_energy_drift_hartree": exact.maximum_energy_drift_hartree,
        "exact_grid_reversal_error": grid_reversal_error,
        "exact_grid_boundary_probability": boundary_probability,
        "grid_refinement_error_dt_004": refinement_errors[0],
        "grid_refinement_error_dt_002": refinement_errors[1],
        "grid_refinement_error_dt_001": refinement_errors[2],
        "grid_order_coarse": grid_order_coarse,
        "grid_order_fine": grid_order_fine,
        "doublet_grid_norm_drift": doublet_exact.maximum_norm_drift,
        "singlet_triplet_grid_norm_drift": singlet_triplet_exact.maximum_norm_drift,
        "overlap_quadrature_error": overlap_quadrature_error,
        "hamiltonian_quadrature_error": hamiltonian_quadrature_error,
        "reduction_overlap_error": reduction_overlap_error,
        "reduction_hamiltonian_error": reduction_hamiltonian_error,
        "reduction_metric_error": reduction_metric_error,
        "reduction_rhs_error": reduction_rhs_error,
        "reduction_velocity_error": reduction_velocity_error,
        "metric_minimum_eigenvalue": metric_minimum_eigenvalue,
        "metric_linear_residual_relative": metric_system.solve_receipt.linear_residual_relative,
        "tdvp_reversal_error": tdvp_reversal_error,
        "gauge_velocity_error": gauge_velocity_error,
        "packet_permutation_velocity_error": permutation_velocity_error,
        "signed_coordinate_permutation_velocity_error": coordinate_velocity_error,
        "one_packet_norm_drift": one_packet.maximum_norm_drift,
        "one_packet_energy_drift_hartree": one_packet.maximum_energy_drift_hartree,
        "one_packet_exact_wavefunction_error": one_packet_exact_error,
        "controlled_exact_wavefunction_error": controlled_exact_error,
        "one_packet_density_error": one_packet_density_error,
        "controlled_density_error": controlled_density_error,
        "adaptive_wavefunction_improvement_ratio": controlled_exact_error / one_packet_exact_error,
        "adaptive_density_improvement_ratio": controlled_density_error / one_packet_density_error,
        "best_candidate_residual_capture": max(item.residual_capture for item in candidate_evaluations),
        "duplicate_candidate_novelty": duplicate_evaluation.novelty,
        "prune_projection_loss": prune_event.projection.relative_projection_loss,
        "prune_energy_jump_hartree": abs(prune_event.projection.energy_jump_hartree),
        "merge_projection_loss": merge_event.projection.relative_projection_loss,
        "merge_energy_jump_hartree": abs(merge_event.projection.energy_jump_hartree),
        "spawn_projection_loss": spawn_event.projection.relative_projection_loss,
        "dormant_shape_drift": dormant_shape_drift,
        "inactive_step_norm_drift": abs(inactive_step.norm_change),
        "no_event_reduction_error": no_event_reduction_error,
        "zero_soc_score_error": zero_soc_score_error,
    }
    checks = {
        "validation_schema_is_v0260": MULTIDIMENSIONAL_VALIDATION_SCHEMA_V260.endswith("v0.26.0"),
        "exact_grid_schema_is_v0260": EXACT_GRID_SCHEMA_V260.endswith("v0.26.0"),
        "two_state_model_is_2d": model.ndim == 2 and model.nstate == 2,
        "soc_opens_origin_gap": abs(origin_gap - 0.02) < 2.0e-13,
        "zero_soc_restores_exact_ci": abs(zero_origin_gap) < 2.0e-13,
        "analytic_derivative_matches_finite_difference": derivative_error < thresholds["model_derivative_error"],
        "coordinate_transformation_is_exact": coordinate_transform_error < thresholds["coordinate_transform_error"],
        "doublet_model_is_complete": doublet_model.nstate == 4 and doublet_model.complete_spin_manifold,
        "doublet_time_reversal_is_preserved": doublet_time_reversal_error < thresholds["time_reversal_symmetry_error"],
        "kramers_pairs_are_degenerate": kramers_splitting < thresholds["kramers_splitting"],
        "singlet_triplet_model_is_complete": singlet_triplet_model.nstate == 5 and singlet_triplet_model.complete_spin_manifold,
        "singlet_projector_rank_is_two": round(np.trace(singlet_triplet_model.projectors["singlet"]).real) == 2,
        "triplet_projector_rank_is_three": round(np.trace(singlet_triplet_model.projectors["triplet"]).real) == 3,
        "exact_grid_norm_is_unitary": exact.maximum_norm_drift < thresholds["grid_norm_drift"],
        "exact_grid_energy_is_bounded": exact.maximum_energy_drift_hartree < thresholds["grid_energy_drift"],
        "exact_grid_is_time_reversible": grid_reversal_error < thresholds["grid_reversal_error"],
        "exact_grid_boundary_is_clear": boundary_probability < thresholds["grid_boundary_probability"],
        "exact_grid_coarse_order_is_second": grid_order_coarse > thresholds["minimum_grid_order"],
        "exact_grid_fine_order_is_second": grid_order_fine > thresholds["minimum_grid_order"],
        "exact_grid_error_decreases_004_to_002": refinement_errors[1] < refinement_errors[0],
        "exact_grid_error_decreases_002_to_001": refinement_errors[2] < refinement_errors[1],
        "doublet_exact_grid_runs": doublet_exact.maximum_norm_drift < thresholds["grid_norm_drift"],
        "singlet_triplet_exact_grid_runs": singlet_triplet_exact.maximum_norm_drift < thresholds["grid_norm_drift"],
        "doublet_projector_populations_resolve_norm": abs(sum(values[-1] for values in doublet_exact.populations.values()) - 1.0) < 2.0e-12,
        "singlet_triplet_populations_resolve_norm": abs(sum(values[-1] for values in singlet_triplet_exact.populations.values()) - 1.0) < 2.0e-12,
        "analytic_overlap_matches_independent_grid": overlap_quadrature_error < thresholds["quadrature_overlap_error"],
        "analytic_hamiltonian_matches_independent_grid": hamiltonian_quadrature_error < thresholds["quadrature_hamiltonian_error"],
        "one_dimensional_overlap_reduces_to_v0252": reduction_overlap_error < thresholds["one_dimensional_reduction_error"],
        "one_dimensional_hamiltonian_reduces_to_v0252": reduction_hamiltonian_error < thresholds["one_dimensional_reduction_error"],
        "one_dimensional_metric_reduces_to_v0252": reduction_metric_error < thresholds["one_dimensional_reduction_error"],
        "one_dimensional_rhs_reduces_to_v0252": reduction_rhs_error < thresholds["one_dimensional_reduction_error"],
        "one_dimensional_velocity_reduces_to_v0252": reduction_velocity_error < thresholds["one_dimensional_reduction_error"],
        "multidimensional_metric_is_psd": metric_minimum_eigenvalue > thresholds["metric_minimum_eigenvalue"],
        "multidimensional_metric_solve_is_accurate": metric_system.solve_receipt.linear_residual_relative < thresholds["metric_linear_residual"],
        "implicit_midpoint_is_reversible": tdvp_reversal_error < thresholds["tdvp_reversal_error"],
        "implicit_midpoint_norm_is_conserved": abs(forward.norm_change) < thresholds["tdvp_norm_drift"],
        "implicit_midpoint_energy_is_conserved": abs(forward.energy_change_hartree) < thresholds["tdvp_energy_drift"],
        "constant_gauge_velocity_is_covariant": gauge_velocity_error < thresholds["covariance_error"],
        "packet_permutation_velocity_is_covariant": permutation_velocity_error < thresholds["covariance_error"],
        "signed_coordinate_permutation_is_covariant": coordinate_velocity_error < thresholds["covariance_error"],
        "one_packet_trajectory_norm_is_conserved": one_packet.maximum_norm_drift < thresholds["tdvp_norm_drift"],
        "one_packet_trajectory_energy_is_conserved": one_packet.maximum_energy_drift_hartree < thresholds["tdvp_energy_drift"],
        "one_packet_matches_exact_short_time": one_packet_exact_error < thresholds["maximum_exact_wavefunction_error"],
        "controlled_basis_matches_exact_short_time": controlled_exact_error < thresholds["maximum_exact_wavefunction_error"],
        "one_packet_density_matches_exact": one_packet_density_error < thresholds["maximum_density_error"],
        "controlled_density_matches_exact": controlled_density_error < thresholds["maximum_density_error"],
        "adaptive_basis_improves_wavefunction_error": controlled_exact_error / one_packet_exact_error < thresholds["adaptive_improvement_ratio"],
        "adaptive_basis_improves_density_error": controlled_density_error / one_packet_density_error < thresholds["adaptive_improvement_ratio"],
        "controlled_trajectory_spawns_once": controlled.event_counts["spawn"] == 1,
        "controlled_trajectory_obeys_packet_cap": controlled.maximum_packet_count <= controlled.settings.maximum_packet_count,
        "candidate_set_covers_all_signed_axes": len(candidates) == 4 * initial.ndim * initial.ngaussian,
        "at_least_one_residual_candidate_is_admitted": any(item.admitted for item in candidate_evaluations),
        "highest_score_candidate_is_selected": spawn_event.selected_candidate.canonical_key() == sorted([item for item in candidate_evaluations if item.admitted], key=lambda item: (-item.residual_capture, item.candidate.canonical_key()))[0].candidate.canonical_key(),
        "duplicate_candidate_is_rejected": duplicate_evaluation.admitted is False,
        "duplicate_candidate_has_zero_novelty": duplicate_evaluation.novelty < 2.0e-12,
        "duplicate_candidate_fires_rank_gate": "rank-deficient-enlarged-basis" in duplicate_evaluation.rejection_reasons,
        "spawn_projection_is_exact": spawn_event.projection.relative_projection_loss < thresholds["projection_loss"],
        "spawn_newborn_coefficient_is_exactly_zero": np.max(np.abs(spawn_event.after.coefficients[-1])) == 0.0,
        "spawn_stable_id_is_monotone": spawn_event.added_packet_id == "g000001",
        "spawn_newborn_age_is_zero": spawn_event.packet_ages_after[-1] == 0,
        "prune_event_is_accepted": prune_event.event_kind == "prune",
        "prune_removes_small_packet": prune_event.removed_packet_id == "g000001",
        "prune_projection_loss_is_bounded": prune_event.projection.relative_projection_loss < thresholds["projection_loss"],
        "prune_energy_jump_is_bounded": abs(prune_event.projection.energy_jump_hartree) < thresholds["projection_energy_jump"],
        "merge_event_is_accepted": merge_event.event_kind == "merge",
        "merge_removes_one_packet": merge_event.after.ngaussian == merge_event.before.ngaussian - 1,
        "merge_projection_loss_is_bounded": merge_event.projection.relative_projection_loss < thresholds["projection_loss"],
        "merge_energy_jump_is_bounded": abs(merge_event.projection.energy_jump_hartree) < thresholds["projection_energy_jump"],
        "newborn_shape_starts_dormant": dormant_mask.tolist() == [True, False],
        "newborn_activation_requires_metric_gate": activated_mask.tolist() == [True, True],
        "dormant_shape_is_exactly_frozen": dormant_shape_drift < thresholds["dormant_shape_drift"],
        "dormant_coefficient_step_conserves_norm": abs(inactive_step.norm_change) < thresholds["tdvp_norm_drift"],
        "no_event_controlled_path_reduces_to_tdvp": no_event_reduction_error < thresholds["no_event_reduction_error"],
        "zero_soc_toggle_preserves_scores": zero_soc_score_error < thresholds["zero_soc_error"],
        "full_width_matrix_claim_remains_false": V260_MULTIDIMENSIONAL_TDVP_CLAIMS["full_correlated_width_matrices_validated"] is False,
        "coordinate_dependent_gauge_claim_remains_false": V260_MULTIDIMENSIONAL_TDVP_CLAIMS["coordinate_dependent_electronic_gauge_covariance_validated"] is False,
        "general_aims_claim_remains_false": V260_MULTIDIMENSIONAL_BASIS_CLAIMS["full_aims_branching_validated"] is False,
        "real_pyscf_trajectory_claim_remains_false": V260_MULTIDIMENSIONAL_BASIS_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
        "absorbing_boundary_claim_remains_false": V260_EXACT_GRID_CLAIMS["absorbing_boundary_conditions_validated"] is False,
        "exact_oracle_independence_claim_is_true": V260_EXACT_GRID_CLAIMS["independent_of_gaussian_tdvp_implementation"] is True,
    }
    evidence = MultidimensionalValidationEvidenceV260(
        metrics=metrics,
        thresholds=thresholds,
        checks={name: bool(value) for name, value in checks.items()},
        exact_grid_fingerprint=exact.fingerprint(),
        controlled_trajectory_fingerprint=controlled.fingerprint(),
    )
    return evidence.validate()


def save_multidimensional_validation_evidence_v260(path):
    evidence = run_multidimensional_validation_evidence_v260()
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(evidence.as_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return evidence
