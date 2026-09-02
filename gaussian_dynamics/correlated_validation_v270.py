"""Independent scientific validation for the v0.27.0 correlated-width release.

The evidence deliberately mixes implementation-level invariants with oracles
that do not use the correlated moment engine: dense FFT quadrature and the
closed Riccati equations for a general Gaussian in a quadratic Hamiltonian.
"""

from dataclasses import dataclass, replace
import hashlib
import itertools
import json

import numpy as np
from scipy.integrate import solve_ivp

from .correlated_basis_adaptation_v270 import (
    CORRELATED_BASIS_SCHEMA_V270,
    ControlledCorrelatedBasisSettingsV270,
    CorrelatedSpawnCandidateV270,
    V270_MULTIDIMENSIONAL_BASIS_CLAIMS,
    adapt_correlated_basis_once_v270,
    evaluate_correlated_spawn_candidate_v270,
    generate_correlated_spawn_candidates_v270,
    metric_compatible_activation_mask_v270,
)
from .correlated_gaussian_tdvp_v270 import (
    CORRELATED_TDVP_SCHEMA_V270,
    CorrelatedGaussianSpinorStateV270,
    V270_CORRELATED_TDVP_CLAIMS,
    build_correlated_gaussian_matrices_v270,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    correlated_moment_table_v270,
    correlated_reduced_density_v270,
    correlated_variational_energy_v270,
    cross_correlated_gaussian_data_v270,
    evaluate_correlated_state_v270,
    exp_frechet_symmetric_v270,
    exp_symmetric_v270,
    gauge_correlated_velocity_v270,
    log_spd_v270,
    pack_correlated_parameters_v270,
    permute_correlated_velocity_v270,
    rotate_correlated_velocity_v270,
    run_correlated_tdvp_v270,
    smat_v270,
    state_from_correlated_parameters_v270,
    symmetric_basis_v270,
    svec_v270,
)
from .multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    build_multidimensional_gaussian_matrices_v260,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
    multidimensional_reduced_density_v260,
    multidimensional_variational_energy_v260,
    pack_multidimensional_parameters_v260,
)
from .multidimensional_soc_v260 import (
    QuadraticSpinHamiltonianNDV260,
    UniformGrid2DV260,
    two_state_ci_soc_model_v260,
)


CORRELATED_VALIDATION_SCHEMA_V270 = "gnd-correlated-width-validation-v0.27.0"


def _sha256_v270(value):
    def canonical(item):
        if isinstance(item, np.generic):
            return canonical(item.item())
        if isinstance(item, complex):
            return [float(item.real), float(item.imag)]
        if isinstance(item, np.ndarray):
            return canonical(item.tolist())
        if isinstance(item, dict):
            return {str(key): canonical(value) for key, value in sorted(item.items())}
        if isinstance(item, (tuple, list)):
            return [canonical(value) for value in item]
        return item

    payload = json.dumps(
        canonical(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _random_unitary_v270(dimension, seed):
    rng = np.random.default_rng(int(seed))
    matrix = rng.normal(size=(dimension, dimension)) + 1.0j * rng.normal(
        size=(dimension, dimension)
    )
    q, r = np.linalg.qr(matrix)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0.0, phases / np.abs(phases), 1.0)
    return q @ np.diag(phases.conj())


def _correlated_state_v270():
    return CorrelatedGaussianSpinorStateV270(
        q=[[-0.25, 0.10]],
        p=[[3.0, 0.20]],
        width_matrices=[[[1.70, 0.25], [0.25, 2.60]]],
        chirp_matrices=[[[0.02, 0.03], [0.03, -0.01]]],
        coefficients=[[1.0, 0.0]],
    ).normalized()


def _quadrature_state_v270():
    return CorrelatedGaussianSpinorStateV270(
        q=[[-0.7, 0.2], [0.8, -0.3]],
        p=[[0.5, -0.2], [-0.3, 0.4]],
        width_matrices=[
            [[1.10, 0.23], [0.23, 0.90]],
            [[0.80, -0.17], [-0.17, 1.20]],
        ],
        chirp_matrices=[
            [[0.10, 0.08], [0.08, -0.04]],
            [[-0.05, 0.06], [0.06, 0.08]],
        ],
        coefficients=[[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()


def _parameter_rotation_matrix_v270(state, orthogonal):
    """Orthogonal packed-coordinate map induced by ``R' = O R``."""

    state = state.validate(require_normalized=False)
    orthogonal = np.asarray(orthogonal, dtype=float)
    size = state.parameter_count
    transform = np.zeros((size, size), dtype=float)
    coefficient_size = state.ngaussian * state.nstate
    transform[: 2 * coefficient_size, : 2 * coefficient_size] = np.eye(
        2 * coefficient_size
    )
    cursor = 2 * coefficient_size
    for family_size in (state.ndim, state.ndim):
        for packet in range(state.ngaussian):
            block = slice(cursor + packet * family_size, cursor + (packet + 1) * family_size)
            transform[block, block] = orthogonal
        cursor += state.ngaussian * family_size
    matrix_transform = np.column_stack(
        [svec_v270(orthogonal @ basis @ orthogonal.T) for basis in symmetric_basis_v270(state.ndim)]
    )
    for _ in range(2):
        for packet in range(state.ngaussian):
            block = slice(
                cursor + packet * state.symmetric_size,
                cursor + (packet + 1) * state.symmetric_size,
            )
            transform[block, block] = matrix_transform
        cursor += state.ngaussian * state.symmetric_size
    return transform


def _independent_fft_quadrature_v270(state, model):
    """Direct grid/FFT S, H, and fourth moments with no analytic moment calls."""

    grid = UniformGrid2DV260.from_bounds((-9.0, 9.0), (-9.0, 9.0), (192, 192))
    points = grid.mesh()
    gaussians = []
    for packet in range(state.ngaussian):
        displacement = points - state.q[packet]
        width = state.width_matrices[packet]
        chirp = state.chirp_matrices[packet]
        normalization = float((np.linalg.det(width) / np.pi**state.ndim) ** 0.25)
        exponent = (
            -0.5 * np.einsum("...a,ab,...b->...", displacement, width, displacement)
            + 0.5j * np.einsum("...a,ab,...b->...", displacement, chirp, displacement)
            + 1.0j * np.einsum("...a,a->...", displacement, state.p[packet])
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
        kinetic = np.fft.ifftn(kinetic_energy * np.fft.fftn(gaussians[packet_j]))
        for electronic_j in range(state.nstate):
            ket = packet_j * state.nstate + electronic_j
            action = potential[..., :, electronic_j] * gaussians[packet_j][..., None]
            action[..., electronic_j] += kinetic
            for packet_i in range(state.ngaussian):
                for electronic_i in range(state.nstate):
                    bra = packet_i * state.nstate + electronic_i
                    if electronic_i == electronic_j:
                        overlap_grid[bra, ket] = (
                            np.vdot(gaussians[packet_i], gaussians[packet_j])
                            * grid.volume_element
                        )
                    hamiltonian_grid[bra, ket] = (
                        np.vdot(gaussians[packet_i], action[..., electronic_i])
                        * grid.volume_element
                    )
    overlap_analytic, hamiltonian_analytic = build_correlated_gaussian_matrices_v270(
        state, model
    )
    pair = (0, 1)
    data = cross_correlated_gaussian_data_v270(
        state.q[pair[0]], state.p[pair[0]], state.width_matrices[pair[0]],
        state.chirp_matrices[pair[0]], state.q[pair[1]], state.p[pair[1]],
        state.width_matrices[pair[1]], state.chirp_matrices[pair[1]],
    )
    overlap, moments = correlated_moment_table_v270(*data, maximum_degree=4)
    product = np.conj(gaussians[pair[0]]) * gaussians[pair[1]]
    moment_errors = []
    for powers in itertools.product(range(5), repeat=2):
        if sum(powers) <= 4:
            direct = np.sum(
                product * points[..., 0] ** powers[0] * points[..., 1] ** powers[1]
            ) * grid.volume_element
            moment_errors.append(abs(direct - overlap * moments[powers]))
    direct_norm = 0.0
    wavefunction = np.einsum("ixy,ia->axy", gaussians, state.coefficients, optimize=True)
    direct_norm = float(np.real(np.vdot(wavefunction, wavefunction)) * grid.volume_element)
    return {
        "overlap_error": float(np.max(np.abs(overlap_grid - overlap_analytic))),
        "hamiltonian_error": float(np.max(np.abs(hamiltonian_grid - hamiltonian_analytic))),
        "moment_error": float(max(moment_errors)),
        "direct_norm_error": abs(direct_norm - state.generalized_norm),
    }


def _riccati_validation_v270():
    angle = 0.43
    rotation = np.asarray(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    H2 = rotation @ np.diag([0.08, 0.23]) @ rotation.T
    mass = np.asarray([[5.0, 0.4], [0.4, 7.0]])
    inverse_mass = np.linalg.inv(mass)
    H1 = np.asarray([[[0.03]], [[-0.02]]])
    model = QuadraticSpinHamiltonianNDV260(
        mass, [[0.0]], H1, H2[:, :, None, None], label="rotated scalar harmonic oracle"
    ).validate()
    width = np.asarray([[2.0, 0.0], [0.0, 2.0]])
    chirp = np.asarray([[0.10, 0.0], [0.0, -0.05]])
    state = CorrelatedGaussianSpinorStateV270(
        [[0.2, -0.3]], [[0.4, -0.2]], [width], [chirp], [[1.0]]
    ).normalized()
    system = build_correlated_metric_system_v270(state, model)
    velocity = system.velocity
    qdot = velocity[2:4]
    pdot = velocity[4:6]
    edot = smat_v270(velocity[6:9], 2)
    gammadot = exp_frechet_symmetric_v270(state.log_width_matrices[0], edot)
    bdot = smat_v270(velocity[9:12], 2)
    precision = width - 1.0j * chirp
    precision_dot = -1.0j * precision @ inverse_mass @ precision + 2.0j * H2
    exact_qdot = inverse_mass @ state.p[0]
    exact_pdot = -(H1[:, 0, 0].real + 2.0 * H2 @ state.q[0])
    exact_gammadot = precision_dot.real
    exact_bdot = -precision_dot.imag

    initial_physical = np.concatenate(
        (state.q[0], state.p[0], width.reshape(-1), chirp.reshape(-1))
    )

    def physical_rhs(_time, values):
        q = values[:2]
        p = values[2:4]
        gamma = values[4:8].reshape(2, 2)
        beta = values[8:12].reshape(2, 2)
        K = gamma - 1.0j * beta
        Kdot = -1.0j * K @ inverse_mass @ K + 2.0j * H2
        return np.concatenate(
            (
                inverse_mass @ p,
                -(H1[:, 0, 0].real + 2.0 * H2 @ q),
                Kdot.real.reshape(-1),
                (-Kdot.imag).reshape(-1),
            )
        )

    reference = solve_ivp(
        physical_rhs, (0.0, 0.2), initial_physical, rtol=1.0e-12, atol=1.0e-14
    ).y[:, -1]
    errors = []
    trajectories = []
    for dt in (0.04, 0.02, 0.01):
        trajectory = run_correlated_tdvp_v270(state, model, dt, round(0.2 / dt))
        final = trajectory.final_state
        physical = np.concatenate(
            (
                final.q[0], final.p[0], final.width_matrices[0].reshape(-1),
                final.chirp_matrices[0].reshape(-1),
            )
        )
        errors.append(float(np.max(np.abs(physical - reference))))
        trajectories.append(trajectory)
    orders = (
        float(np.log2(errors[0] / errors[1])),
        float(np.log2(errors[1] / errors[2])),
    )
    return {
        "q_velocity_error": float(np.max(np.abs(qdot - exact_qdot))),
        "p_velocity_error": float(np.max(np.abs(pdot - exact_pdot))),
        "width_velocity_error": float(np.max(np.abs(gammadot - exact_gammadot))),
        "chirp_velocity_error": float(np.max(np.abs(bdot - exact_bdot))),
        "analytic_correlation_rate": float(
            max(abs(exact_gammadot[0, 1]), abs(exact_bdot[0, 1]))
        ),
        "errors": tuple(errors),
        "orders": orders,
        "maximum_norm_drift": max(item.maximum_norm_drift for item in trajectories),
        "maximum_energy_drift": max(
            item.maximum_energy_drift_hartree for item in trajectories
        ),
        "metric_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(system.metric))),
        "metric_residual": system.solve_receipt.linear_residual_relative,
        "trajectory_fingerprint": trajectories[-1].fingerprint(),
    }


@dataclass(frozen=True)
class CorrelatedValidationEvidenceV270:
    metrics: dict
    thresholds: dict
    checks: dict
    trajectory_fingerprint: str
    lifecycle_fingerprint: str

    @property
    def passed(self):
        return all(self.checks.values())

    @property
    def check_count(self):
        return len(self.checks)

    def validate(self):
        if self.check_count != 100:
            raise ValueError(
                f"v0.27.0 validation evidence requires exactly 100 gates, found {self.check_count}."
            )
        if any(type(value) is not bool for value in self.checks.values()):
            raise TypeError("every v0.27.0 validation gate must be a native Boolean.")
        if not self.passed:
            failed = ", ".join(name for name, value in self.checks.items() if not value)
            raise ValueError("v0.27.0 validation failed: " + failed)
        for mapping in (self.metrics, self.thresholds):
            for name, value in mapping.items():
                if not np.isfinite(float(value)):
                    raise ValueError(f"validation scalar {name!r} is non-finite.")
        for fingerprint in (self.trajectory_fingerprint, self.lifecycle_fingerprint):
            if len(fingerprint) != 64 or any(character not in "0123456789abcdef" for character in fingerprint):
                raise ValueError("validation evidence fingerprints must be SHA-256 strings.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": CORRELATED_VALIDATION_SCHEMA_V270,
            "metrics": {name: float(value) for name, value in self.metrics.items()},
            "thresholds": {name: float(value) for name, value in self.thresholds.items()},
            "checks": dict(self.checks),
            "passed": self.passed,
            "check_count": self.check_count,
            "fingerprints": {
                "riccati_trajectory": self.trajectory_fingerprint,
                "lifecycle": self.lifecycle_fingerprint,
            },
            "claims": {
                "tdvp": dict(V270_CORRELATED_TDVP_CLAIMS),
                "basis": dict(V270_MULTIDIMENSIONAL_BASIS_CLAIMS),
            },
        }

    def fingerprint(self):
        return _sha256_v270(self.as_dict())


def run_correlated_validation_evidence_v270():
    thresholds = {
        "algebra_error": 3.0e-10,
        "frechet_error": 2.0e-8,
        "quadrature_overlap_error": 3.0e-10,
        "quadrature_hamiltonian_error": 3.0e-9,
        "quadrature_moment_error": 3.0e-10,
        "reduction_error": 3.0e-11,
        "instantaneous_oracle_error": 3.0e-10,
        "maximum_riccati_error": 1.0e-5,
        "minimum_midpoint_order": 1.95,
        "metric_minimum_eigenvalue": -3.0e-10,
        "metric_residual": 3.0e-9,
        "norm_drift": 1.0e-8,
        "energy_drift": 1.0e-7,
        "covariance_error": 3.0e-8,
        "step_covariance_error": 3.0e-7,
        "time_reversal_error": 3.0e-8,
        "projection_loss": 2.0e-7,
        "projection_energy_jump": 2.0e-6,
        "dormant_shape_drift": 2.0e-13,
        "candidate_score_covariance": 3.0e-10,
    }

    # Symmetric-matrix coordinates, positivity, and Frechet calculus.
    symmetric = np.asarray([[0.37, -0.22], [-0.22, 0.81]])
    other = np.asarray([[-0.14, 0.31], [0.31, 0.26]])
    roundtrip_error = float(np.max(np.abs(smat_v270(svec_v270(symmetric), 2) - symmetric)))
    frobenius_error = abs(float(np.vdot(svec_v270(symmetric), svec_v270(other))) - float(np.vdot(symmetric, other)))
    basis = symmetric_basis_v270(2)
    basis_gram = np.asarray([[np.vdot(left, right) for right in basis] for left in basis])
    basis_error = float(np.max(np.abs(basis_gram - np.eye(3))))
    angle = 0.371
    rotation = np.asarray([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    reflection = np.asarray([[1.0, 0.0], [0.0, -1.0]])
    matrix_rotation = np.column_stack([svec_v270(rotation @ item @ rotation.T) for item in basis])
    matrix_reflection = np.column_stack([svec_v270(reflection @ item @ reflection.T) for item in basis])
    svec_rotation_orthogonality = float(np.max(np.abs(matrix_rotation.T @ matrix_rotation - np.eye(3))))
    svec_reflection_orthogonality = float(np.max(np.abs(matrix_reflection.T @ matrix_reflection - np.eye(3))))
    log_matrix = np.asarray([[0.31, -0.18], [-0.18, -0.27]])
    width = exp_symmetric_v270(log_matrix)
    log_exp_error = float(np.max(np.abs(log_spd_v270(width) - log_matrix)))
    exp_log_error = float(np.max(np.abs(exp_symmetric_v270(log_spd_v270(width)) - width)))
    direction = np.asarray([[0.23, -0.11], [-0.11, 0.36]])
    epsilon = 2.0e-6
    frechet = exp_frechet_symmetric_v270(log_matrix, direction)
    finite_difference = (
        exp_symmetric_v270(log_matrix + epsilon * direction)
        - exp_symmetric_v270(log_matrix - epsilon * direction)
    ) / (2.0 * epsilon)
    frechet_error = float(np.max(np.abs(frechet - finite_difference)))
    frechet_symmetry_error = float(np.max(np.abs(frechet - frechet.T)))
    rotated_frechet = exp_frechet_symmetric_v270(
        rotation @ log_matrix @ rotation.T, rotation @ direction @ rotation.T
    )
    frechet_rotation_error = float(
        np.max(np.abs(rotated_frechet - rotation @ frechet @ rotation.T))
    )
    extreme_width = exp_symmetric_v270(np.asarray([[-16.0, 0.0], [0.0, 12.0]]))
    extreme_minimum_eigenvalue = float(np.min(np.linalg.eigvalsh(extreme_width)))
    log_rotation_error = float(
        np.max(
            np.abs(
                log_spd_v270(rotation @ width @ rotation.T)
                - rotation @ log_spd_v270(width) @ rotation.T
            )
        )
    )

    # Independent moments and matrix elements.
    state = _correlated_state_v270()
    quadrature_state = _quadrature_state_v270()
    model = two_state_ci_soc_model_v260()
    quadrature_model = two_state_ci_soc_model_v260(
        mass_au=(900.0, 700.0), kappa=0.003, coupling=0.004,
        frequencies=(0.02, 0.025), soc_scale=0.002,
    )
    self_data = cross_correlated_gaussian_data_v270(
        state.q[0], state.p[0], state.width_matrices[0], state.chirp_matrices[0],
        state.q[0], state.p[0], state.width_matrices[0], state.chirp_matrices[0],
    )
    cross_01 = cross_correlated_gaussian_data_v270(
        quadrature_state.q[0], quadrature_state.p[0], quadrature_state.width_matrices[0],
        quadrature_state.chirp_matrices[0], quadrature_state.q[1], quadrature_state.p[1],
        quadrature_state.width_matrices[1], quadrature_state.chirp_matrices[1],
    )
    cross_10 = cross_correlated_gaussian_data_v270(
        quadrature_state.q[1], quadrature_state.p[1], quadrature_state.width_matrices[1],
        quadrature_state.chirp_matrices[1], quadrature_state.q[0], quadrature_state.p[0],
        quadrature_state.width_matrices[0], quadrature_state.chirp_matrices[0],
    )
    overlap_self_error = abs(self_data[0] - 1.0)
    overlap_conjugation_error = abs(cross_01[0] - np.conj(cross_10[0]))
    rotated_quadrature_state = quadrature_state.coordinate_rotated(rotation)
    rotated_cross = cross_correlated_gaussian_data_v270(
        rotated_quadrature_state.q[0], rotated_quadrature_state.p[0],
        rotated_quadrature_state.width_matrices[0], rotated_quadrature_state.chirp_matrices[0],
        rotated_quadrature_state.q[1], rotated_quadrature_state.p[1],
        rotated_quadrature_state.width_matrices[1], rotated_quadrature_state.chirp_matrices[1],
    )
    reflected_quadrature_state = quadrature_state.coordinate_rotated(reflection)
    reflected_cross = cross_correlated_gaussian_data_v270(
        reflected_quadrature_state.q[0], reflected_quadrature_state.p[0],
        reflected_quadrature_state.width_matrices[0], reflected_quadrature_state.chirp_matrices[0],
        reflected_quadrature_state.q[1], reflected_quadrature_state.p[1],
        reflected_quadrature_state.width_matrices[1], reflected_quadrature_state.chirp_matrices[1],
    )
    overlap_rotation_error = abs(rotated_cross[0] - cross_01[0])
    overlap_reflection_error = abs(reflected_cross[0] - cross_01[0])
    covariance_symmetry_error = float(np.max(np.abs(cross_01[2] - cross_01[2].T)))
    moment_zero_error = abs(correlated_moment_table_v270(*cross_01, maximum_degree=4)[1][(0, 0)] - 1.0)
    quadrature = _independent_fft_quadrature_v270(quadrature_state, quadrature_model)
    overlap_matrix, hamiltonian_matrix = build_correlated_gaussian_matrices_v270(
        quadrature_state, quadrature_model
    )
    overlap_hermiticity = float(np.max(np.abs(overlap_matrix - overlap_matrix.conj().T)))
    hamiltonian_hermiticity = float(np.max(np.abs(hamiltonian_matrix - hamiltonian_matrix.conj().T)))
    energy = correlated_variational_energy_v270(quadrature_state, quadrature_model)
    density = correlated_reduced_density_v270(quadrature_state)
    density_hermiticity = float(np.max(np.abs(density - density.conj().T)))
    density_trace_error = abs(float(np.trace(density).real) - 1.0)

    # Exact reductions to the diagonal v0.26.0 manifold and its 1D special case.
    diagonal_state = DiagonalGaussianSpinorStateV260(
        q=[[-0.7, 0.2], [0.8, -0.3]], p=[[0.5, -0.2], [-0.3, 0.4]],
        widths=[[1.1, 0.9], [0.8, 1.2]], chirps=[[0.1, -0.04], [-0.05, 0.08]],
        coefficients=[[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    embedded = CorrelatedGaussianSpinorStateV270.from_diagonal_v260(diagonal_state)
    recovered = embedded.to_diagonal_v260()
    conversion_error = float(
        np.max(
            np.abs(
                pack_multidimensional_parameters_v260(recovered)
                - pack_multidimensional_parameters_v260(diagonal_state)
            )
        )
    )
    old_matrices_2d = build_multidimensional_gaussian_matrices_v260(diagonal_state, quadrature_model)
    new_matrices_2d = build_correlated_gaussian_matrices_v270(embedded, quadrature_model)
    reduction_2d_overlap = float(np.max(np.abs(old_matrices_2d[0] - new_matrices_2d[0])))
    reduction_2d_hamiltonian = float(np.max(np.abs(old_matrices_2d[1] - new_matrices_2d[1])))
    reduction_2d_energy = abs(
        multidimensional_variational_energy_v260(diagonal_state, quadrature_model)
        - correlated_variational_energy_v270(embedded, quadrature_model)
    )
    reduction_2d_density = float(
        np.max(
            np.abs(
                multidimensional_reduced_density_v260(diagonal_state)
                - correlated_reduced_density_v270(embedded)
            )
        )
    )
    H0_1d = np.asarray([[0.01, 0.002j], [-0.002j, -0.01]])
    H1_1d = np.asarray([[[0.003, 0.004], [0.004, -0.002]]])
    H2_1d = np.asarray([[[[0.0005, 0.0], [0.0, 0.0007]]]])
    model_1d = QuadraticSpinHamiltonianNDV260([[900.0]], H0_1d, H1_1d, H2_1d).validate()
    state_1d_old = DiagonalGaussianSpinorStateV260(
        [[-0.7], [0.8]], [[0.5], [-0.3]], [[1.1], [0.8]], [[0.1], [-0.05]],
        [[0.7, 0.2j], [0.15 - 0.1j, 0.3]],
    ).normalized()
    state_1d_new = CorrelatedGaussianSpinorStateV270.from_diagonal_v260(state_1d_old)
    old_matrices_1d = build_multidimensional_gaussian_matrices_v260(state_1d_old, model_1d)
    new_matrices_1d = build_correlated_gaussian_matrices_v270(state_1d_new, model_1d)
    old_metric_1d = build_multidimensional_metric_system_v260(state_1d_old, model_1d)
    new_metric_1d = build_correlated_metric_system_v270(state_1d_new, model_1d)
    reduction_1d_overlap = float(np.max(np.abs(old_matrices_1d[0] - new_matrices_1d[0])))
    reduction_1d_hamiltonian = float(np.max(np.abs(old_matrices_1d[1] - new_matrices_1d[1])))
    reduction_1d_metric = float(np.max(np.abs(old_metric_1d.metric - new_metric_1d.metric)))
    reduction_1d_rhs = float(np.max(np.abs(old_metric_1d.rhs - new_metric_1d.rhs)))
    reduction_1d_velocity = float(np.max(np.abs(old_metric_1d.velocity - new_metric_1d.velocity)))
    old_step_1d = multidimensional_implicit_midpoint_step_v260(state_1d_old, model_1d, 0.001)
    new_step_1d = correlated_implicit_midpoint_step_v270(state_1d_new, model_1d, 0.001)
    reduction_1d_step = float(
        np.max(
            np.abs(
                pack_multidimensional_parameters_v260(old_step_1d.end)
                - pack_correlated_parameters_v270(new_step_1d.end)
            )
        )
    )

    # Closed Riccati oracle and second-order implicit-midpoint convergence.
    riccati = _riccati_validation_v270()

    # Proper rotations, reflections, electronic gauge, and packet permutations.
    coordinate = np.asarray([-0.13, 0.44])
    rotated_state = state.coordinate_rotated(rotation)
    rotated_model = model.coordinate_rotated(rotation)
    model_rotation_error = float(
        np.max(
            np.abs(
                rotated_model.hamiltonian(coordinate)
                - model.hamiltonian(coordinate @ rotation)
            )
        )
    )
    wavefunction_rotation_error = float(
        np.max(
            np.abs(
                evaluate_correlated_state_v270(rotated_state, coordinate)
                - evaluate_correlated_state_v270(state, coordinate @ rotation)
            )
        )
    )
    base_matrices = build_correlated_gaussian_matrices_v270(state, model)
    rotated_matrices = build_correlated_gaussian_matrices_v270(rotated_state, rotated_model)
    rotation_overlap_error = float(np.max(np.abs(base_matrices[0] - rotated_matrices[0])))
    rotation_hamiltonian_error = float(np.max(np.abs(base_matrices[1] - rotated_matrices[1])))
    rotation_energy_error = abs(
        correlated_variational_energy_v270(state, model)
        - correlated_variational_energy_v270(rotated_state, rotated_model)
    )
    rotation_density_error = float(
        np.max(
            np.abs(
                correlated_reduced_density_v270(state)
                - correlated_reduced_density_v270(rotated_state)
            )
        )
    )
    base_system = build_correlated_metric_system_v270(state, model)
    rotated_system = build_correlated_metric_system_v270(rotated_state, rotated_model)
    packed_rotation = _parameter_rotation_matrix_v270(state, rotation)
    metric_rotation_error = float(
        np.max(np.abs(rotated_system.metric - packed_rotation @ base_system.metric @ packed_rotation.T))
    )
    rhs_rotation_error = float(
        np.max(np.abs(rotated_system.rhs - packed_rotation @ base_system.rhs))
    )
    velocity_rotation_error = float(
        np.max(
            np.abs(
                rotated_system.velocity
                - rotate_correlated_velocity_v270(state, base_system.velocity, rotation)
            )
        )
    )
    base_step = correlated_implicit_midpoint_step_v270(state, model, 0.002)
    rotated_step = correlated_implicit_midpoint_step_v270(rotated_state, rotated_model, 0.002)
    step_rotation_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(base_step.end.coordinate_rotated(rotation))
                - pack_correlated_parameters_v270(rotated_step.end)
            )
        )
    )
    reflected_state = state.coordinate_rotated(reflection)
    reflected_model = model.coordinate_rotated(reflection)
    model_reflection_error = float(
        np.max(
            np.abs(
                reflected_model.hamiltonian(coordinate)
                - model.hamiltonian(coordinate @ reflection)
            )
        )
    )
    wavefunction_reflection_error = float(
        np.max(
            np.abs(
                evaluate_correlated_state_v270(reflected_state, coordinate)
                - evaluate_correlated_state_v270(state, coordinate @ reflection)
            )
        )
    )
    reflected_matrices = build_correlated_gaussian_matrices_v270(reflected_state, reflected_model)
    reflection_overlap_error = float(np.max(np.abs(base_matrices[0] - reflected_matrices[0])))
    reflection_hamiltonian_error = float(np.max(np.abs(base_matrices[1] - reflected_matrices[1])))
    reflected_system = build_correlated_metric_system_v270(reflected_state, reflected_model)
    velocity_reflection_error = float(
        np.max(
            np.abs(
                reflected_system.velocity
                - rotate_correlated_velocity_v270(state, base_system.velocity, reflection)
            )
        )
    )
    reflected_step = correlated_implicit_midpoint_step_v270(state.coordinate_rotated(reflection), reflected_model, 0.002)
    step_reflection_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(base_step.end.coordinate_rotated(reflection))
                - pack_correlated_parameters_v270(reflected_step.end)
            )
        )
    )
    unitary = _random_unitary_v270(2, 27001)
    gauge_state = quadrature_state.gauge_transformed(unitary)
    gauge_model = quadrature_model.gauge_transformed(unitary)
    gauge_base_system = build_correlated_metric_system_v270(quadrature_state, quadrature_model)
    gauge_system = build_correlated_metric_system_v270(gauge_state, gauge_model)
    gauge_velocity_error = float(
        np.max(
            np.abs(
                gauge_system.velocity
                - gauge_correlated_velocity_v270(
                    quadrature_state, gauge_base_system.velocity, unitary
                )
            )
        )
    )
    gauge_base_step = correlated_implicit_midpoint_step_v270(
        quadrature_state, quadrature_model, 0.001
    )
    gauge_step = correlated_implicit_midpoint_step_v270(gauge_state, gauge_model, 0.001)
    gauge_step_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(gauge_base_step.end.gauge_transformed(unitary))
                - pack_correlated_parameters_v270(gauge_step.end)
            )
        )
    )
    permutation = np.asarray([1, 0])
    permuted_state = quadrature_state.permuted(permutation)
    permutation_system = build_correlated_metric_system_v270(permuted_state, quadrature_model)
    permutation_velocity_error = float(
        np.max(
            np.abs(
                permutation_system.velocity
                - permute_correlated_velocity_v270(
                    quadrature_state, gauge_base_system.velocity, permutation
                )
            )
        )
    )
    permutation_step = correlated_implicit_midpoint_step_v270(
        permuted_state, quadrature_model, 0.001
    )
    permutation_step_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(gauge_base_step.end.permuted(permutation))
                - pack_correlated_parameters_v270(permutation_step.end)
            )
        )
    )
    backward = correlated_implicit_midpoint_step_v270(base_step.end, model, -0.002)
    time_reversal_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(backward.end)
                - pack_correlated_parameters_v270(state)
            )
        )
    )

    # Controlled correlated basis lifecycle and intrinsic-axis covariance.
    candidates = generate_correlated_spawn_candidates_v270(state)
    evaluations = [
        evaluate_correlated_spawn_candidate_v270(state, model, candidate)
        for candidate in candidates
    ]
    direction_errors = []
    eigenvalues, eigenvectors = np.linalg.eigh(state.width_matrices[0])
    for candidate in candidates:
        displacement = (
            candidate.q - state.q[0]
            if candidate.displacement_kind == "position"
            else candidate.p - state.p[0]
        )
        direction = eigenvectors[:, candidate.coordinate_axis]
        direction_errors.append(
            np.linalg.norm(displacement - direction * np.dot(direction, displacement))
        )
    principal_axis_error = float(max(direction_errors))
    spawn_event = adapt_correlated_basis_once_v270(state, model)
    isotropic = CorrelatedGaussianSpinorStateV270(
        state.q, state.p, [2.0 * np.eye(2)], state.chirp_matrices, state.coefficients
    ).normalized()
    near_degenerate = CorrelatedGaussianSpinorStateV270(
        state.q, state.p, [np.diag([2.0, 2.0 + 1.0e-10])],
        state.chirp_matrices, state.coefficients,
    ).normalized()
    rotated_candidates = generate_correlated_spawn_candidates_v270(rotated_state)
    rotated_evaluations = [
        evaluate_correlated_spawn_candidate_v270(rotated_state, rotated_model, candidate)
        for candidate in rotated_candidates
    ]
    candidate_score_rotation_error = float(
        np.max(
            np.abs(
                np.sort([item.residual_capture for item in evaluations])
                - np.sort([item.residual_capture for item in rotated_evaluations])
            )
        )
    )
    rotated_spawn_event = adapt_correlated_basis_once_v270(rotated_state, rotated_model)
    spawn_event_rotation_error = float(
        np.max(
            np.abs(
                pack_correlated_parameters_v270(spawn_event.after.coordinate_rotated(rotation))
                - pack_correlated_parameters_v270(rotated_spawn_event.after)
            )
        )
    )
    duplicate = CorrelatedSpawnCandidateV270(
        state.q[0], state.p[0], state.width_matrices[0], state.chirp_matrices[0],
        0, "position", 0, 1,
    )
    duplicate_evaluation = evaluate_correlated_spawn_candidate_v270(state, model, duplicate)
    zero_soc_model = two_state_ci_soc_model_v260(soc_scale=0.0)
    lifecycle_width = np.asarray([[2.0, 0.25], [0.25, 1.4]])
    lifecycle_chirp = np.asarray([[0.02, 0.03], [0.03, -0.01]])
    prune_state = CorrelatedGaussianSpinorStateV270(
        [[-2.0, 0.0], [2.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [lifecycle_width, lifecycle_width], [lifecycle_chirp, lifecycle_chirp],
        [[1.0, 0.0], [1.0e-7, 0.0]],
    ).normalized()
    prune_event = adapt_correlated_basis_once_v270(
        prune_state, zero_soc_model, packet_ids=("g000000", "g000001"),
        packet_ages=(64, 64), next_packet_serial=2,
    )
    merge_state = CorrelatedGaussianSpinorStateV270(
        [[0.0, 0.0], [1.0e-4, 0.0]], [[0.0, 0.0], [0.0, 0.0]],
        [lifecycle_width, lifecycle_width], [lifecycle_chirp, lifecycle_chirp],
        [[0.7, 0.0], [0.3, 0.0]],
    ).normalized()
    merge_event = adapt_correlated_basis_once_v270(
        merge_state, zero_soc_model, packet_ids=("g000000", "g000001"),
        packet_ages=(2, 2), next_packet_serial=2,
    )
    dormant_state = spawn_event.after
    dormant_mask = metric_compatible_activation_mask_v270(
        dormant_state, model, locked_active_mask=[True, False]
    )
    activated_coefficients = dormant_state.coefficients.copy()
    activated_coefficients[-1, 1] = 0.03
    activated_state = replace(dormant_state, coefficients=activated_coefficients).normalized()
    activated_mask = metric_compatible_activation_mask_v270(
        activated_state, model, locked_active_mask=[True, False]
    )
    inactive_step = correlated_implicit_midpoint_step_v270(
        dormant_state, model, 0.001, active_shape_mask=dormant_mask
    )
    dormant_shape_drift = max(
        float(np.max(np.abs(inactive_step.end.q[-1] - inactive_step.start.q[-1]))),
        float(np.max(np.abs(inactive_step.end.p[-1] - inactive_step.start.p[-1]))),
        float(np.max(np.abs(inactive_step.end.width_matrices[-1] - inactive_step.start.width_matrices[-1]))),
        float(np.max(np.abs(inactive_step.end.chirp_matrices[-1] - inactive_step.start.chirp_matrices[-1]))),
    )

    metrics = {
        "svec_roundtrip_error": roundtrip_error,
        "svec_frobenius_error": frobenius_error,
        "symmetric_basis_error": basis_error,
        "svec_rotation_orthogonality_error": svec_rotation_orthogonality,
        "svec_reflection_orthogonality_error": svec_reflection_orthogonality,
        "log_exp_roundtrip_error": log_exp_error,
        "exp_log_roundtrip_error": exp_log_error,
        "exp_frechet_finite_difference_error": frechet_error,
        "exp_frechet_symmetry_error": frechet_symmetry_error,
        "exp_frechet_rotation_error": frechet_rotation_error,
        "extreme_width_minimum_eigenvalue": extreme_minimum_eigenvalue,
        "log_width_rotation_error": log_rotation_error,
        "self_overlap_error": overlap_self_error,
        "overlap_conjugation_error": overlap_conjugation_error,
        "cross_overlap_magnitude": abs(cross_01[0]),
        "overlap_rotation_error": overlap_rotation_error,
        "overlap_reflection_error": overlap_reflection_error,
        "complex_covariance_symmetry_error": covariance_symmetry_error,
        "zeroth_moment_error": moment_zero_error,
        "fourth_moment_quadrature_error": quadrature["moment_error"],
        "overlap_fft_quadrature_error": quadrature["overlap_error"],
        "hamiltonian_fft_quadrature_error": quadrature["hamiltonian_error"],
        "grid_norm_error": quadrature["direct_norm_error"],
        "overlap_hermiticity_error": overlap_hermiticity,
        "hamiltonian_hermiticity_error": hamiltonian_hermiticity,
        "variational_energy_hartree": energy,
        "density_hermiticity_error": density_hermiticity,
        "density_trace_error": density_trace_error,
        "diagonal_conversion_error": conversion_error,
        "reduction_2d_overlap_error": reduction_2d_overlap,
        "reduction_2d_hamiltonian_error": reduction_2d_hamiltonian,
        "reduction_2d_energy_error": reduction_2d_energy,
        "reduction_2d_density_error": reduction_2d_density,
        "reduction_1d_overlap_error": reduction_1d_overlap,
        "reduction_1d_hamiltonian_error": reduction_1d_hamiltonian,
        "reduction_1d_metric_error": reduction_1d_metric,
        "reduction_1d_rhs_error": reduction_1d_rhs,
        "reduction_1d_velocity_error": reduction_1d_velocity,
        "reduction_1d_step_error": reduction_1d_step,
        "riccati_q_velocity_error": riccati["q_velocity_error"],
        "riccati_p_velocity_error": riccati["p_velocity_error"],
        "riccati_width_velocity_error": riccati["width_velocity_error"],
        "riccati_chirp_velocity_error": riccati["chirp_velocity_error"],
        "riccati_analytic_correlation_rate": riccati["analytic_correlation_rate"],
        "riccati_error_dt_004": riccati["errors"][0],
        "riccati_error_dt_002": riccati["errors"][1],
        "riccati_error_dt_001": riccati["errors"][2],
        "riccati_order_coarse": riccati["orders"][0],
        "riccati_order_fine": riccati["orders"][1],
        "riccati_trajectory_norm_drift": riccati["maximum_norm_drift"],
        "riccati_trajectory_energy_drift": riccati["maximum_energy_drift"],
        "riccati_metric_minimum_eigenvalue": riccati["metric_minimum_eigenvalue"],
        "riccati_metric_residual": riccati["metric_residual"],
        "model_rotation_error": model_rotation_error,
        "wavefunction_rotation_error": wavefunction_rotation_error,
        "matrix_rotation_overlap_error": rotation_overlap_error,
        "matrix_rotation_hamiltonian_error": rotation_hamiltonian_error,
        "rotation_energy_error": rotation_energy_error,
        "rotation_density_error": rotation_density_error,
        "metric_rotation_error": metric_rotation_error,
        "rhs_rotation_error": rhs_rotation_error,
        "velocity_rotation_error": velocity_rotation_error,
        "step_rotation_error": step_rotation_error,
        "model_reflection_error": model_reflection_error,
        "wavefunction_reflection_error": wavefunction_reflection_error,
        "matrix_reflection_overlap_error": reflection_overlap_error,
        "matrix_reflection_hamiltonian_error": reflection_hamiltonian_error,
        "velocity_reflection_error": velocity_reflection_error,
        "step_reflection_error": step_reflection_error,
        "gauge_velocity_error": gauge_velocity_error,
        "gauge_step_error": gauge_step_error,
        "packet_permutation_velocity_error": permutation_velocity_error,
        "packet_permutation_step_error": permutation_step_error,
        "time_reversal_error": time_reversal_error,
        "principal_axis_alignment_error": principal_axis_error,
        "candidate_score_rotation_error": candidate_score_rotation_error,
        "spawn_event_rotation_error": spawn_event_rotation_error,
        "duplicate_candidate_novelty": duplicate_evaluation.novelty,
        "spawn_projection_loss": spawn_event.projection.relative_projection_loss,
        "prune_projection_loss": prune_event.projection.relative_projection_loss,
        "prune_energy_jump": abs(prune_event.projection.energy_jump_hartree),
        "merge_projection_loss": merge_event.projection.relative_projection_loss,
        "merge_energy_jump": abs(merge_event.projection.energy_jump_hartree),
        "dormant_shape_drift": dormant_shape_drift,
    }

    checks = {
        # 1--8: release schema and honest claim boundary.
        "validation_schema_is_v0270": CORRELATED_VALIDATION_SCHEMA_V270.endswith("v0.27.0"),
        "tdvp_schema_is_v0270": CORRELATED_TDVP_SCHEMA_V270.endswith("v0.27.0"),
        "basis_schema_is_v0270": CORRELATED_BASIS_SCHEMA_V270.endswith("v0.27.0"),
        "full_correlated_width_claim_is_true": V270_CORRELATED_TDVP_CLAIMS["full_correlated_width_matrices_validated"] is True,
        "orthogonal_covariance_claim_is_true": V270_CORRELATED_TDVP_CLAIMS["arbitrary_orthogonal_coordinate_covariance_validated"] is True,
        "log_euclidean_claim_is_true": V270_CORRELATED_TDVP_CLAIMS["log_euclidean_positive_width_parameterization_validated"] is True,
        "live_molecular_soc_claim_remains_false": V270_CORRELATED_TDVP_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False,
        "full_aims_claim_remains_false": V270_MULTIDIMENSIONAL_BASIS_CLAIMS["full_aims_branching_validated"] is False,
        # 9--20: full symmetric-matrix coordinate algebra.
        "svec_roundtrip_is_exact": roundtrip_error < thresholds["algebra_error"],
        "svec_is_frobenius_isometric": frobenius_error < thresholds["algebra_error"],
        "symmetric_basis_is_orthonormal": basis_error < thresholds["algebra_error"],
        "proper_rotation_is_orthogonal_in_svec": svec_rotation_orthogonality < thresholds["algebra_error"],
        "reflection_is_orthogonal_in_svec": svec_reflection_orthogonality < thresholds["algebra_error"],
        "log_after_exp_roundtrip_is_exact": log_exp_error < thresholds["algebra_error"],
        "exp_after_log_roundtrip_is_exact": exp_log_error < thresholds["algebra_error"],
        "exp_frechet_matches_finite_difference": frechet_error < thresholds["frechet_error"],
        "exp_frechet_is_symmetric": frechet_symmetry_error < thresholds["algebra_error"],
        "exp_frechet_is_rotation_covariant": frechet_rotation_error < thresholds["algebra_error"],
        "extreme_log_width_remains_spd": extreme_minimum_eigenvalue > 0.0,
        "log_width_commutes_with_rotation": log_rotation_error < thresholds["algebra_error"],
        # 21--35: correlated moments and independent FFT quadrature.
        "normalized_packet_self_overlap_is_one": overlap_self_error < thresholds["algebra_error"],
        "cross_overlap_obeys_conjugation": overlap_conjugation_error < thresholds["algebra_error"],
        "cross_overlap_respects_cauchy_schwarz": abs(cross_01[0]) <= 1.0 + thresholds["algebra_error"],
        "cross_overlap_is_rotation_invariant": overlap_rotation_error < thresholds["algebra_error"],
        "cross_overlap_is_reflection_invariant": overlap_reflection_error < thresholds["algebra_error"],
        "complex_covariance_is_symmetric": covariance_symmetry_error < thresholds["algebra_error"],
        "zeroth_moment_is_one": moment_zero_error < thresholds["algebra_error"],
        "fourth_moments_match_direct_quadrature": quadrature["moment_error"] < thresholds["quadrature_moment_error"],
        "analytic_overlap_matches_fft_quadrature": quadrature["overlap_error"] < thresholds["quadrature_overlap_error"],
        "analytic_hamiltonian_matches_fft_quadrature": quadrature["hamiltonian_error"] < thresholds["quadrature_hamiltonian_error"],
        "grid_wavefunction_norm_matches_overlap_norm": quadrature["direct_norm_error"] < thresholds["quadrature_overlap_error"],
        "analytic_overlap_matrix_is_hermitian": overlap_hermiticity < thresholds["algebra_error"],
        "analytic_hamiltonian_matrix_is_hermitian": hamiltonian_hermiticity < thresholds["algebra_error"],
        "variational_energy_is_finite": np.isfinite(energy),
        "reduced_density_is_hermitian_and_normalized": max(density_hermiticity, density_trace_error) < thresholds["algebra_error"],
        # 36--46: exact inherited-manifold reductions.
        "diagonal_state_conversion_roundtrips": conversion_error < thresholds["reduction_error"],
        "two_dimensional_overlap_reduces_to_v0260": reduction_2d_overlap < thresholds["reduction_error"],
        "two_dimensional_hamiltonian_reduces_to_v0260": reduction_2d_hamiltonian < thresholds["reduction_error"],
        "two_dimensional_energy_reduces_to_v0260": reduction_2d_energy < thresholds["reduction_error"],
        "two_dimensional_density_reduces_to_v0260": reduction_2d_density < thresholds["reduction_error"],
        "one_dimensional_overlap_reduces_to_v0260": reduction_1d_overlap < thresholds["reduction_error"],
        "one_dimensional_hamiltonian_reduces_to_v0260": reduction_1d_hamiltonian < thresholds["reduction_error"],
        "one_dimensional_metric_reduces_to_v0260": reduction_1d_metric < thresholds["reduction_error"],
        "one_dimensional_rhs_reduces_to_v0260": reduction_1d_rhs < thresholds["reduction_error"],
        "one_dimensional_velocity_reduces_to_v0260": reduction_1d_velocity < thresholds["reduction_error"],
        "one_dimensional_midpoint_step_reduces_to_v0260": reduction_1d_step < thresholds["reduction_error"],
        # 47--60: exact quadratic-Gaussian Riccati oracle.
        "centroid_velocity_matches_hamilton_equation": riccati["q_velocity_error"] < thresholds["instantaneous_oracle_error"],
        "momentum_velocity_matches_hamilton_equation": riccati["p_velocity_error"] < thresholds["instantaneous_oracle_error"],
        "full_width_velocity_matches_riccati_equation": riccati["width_velocity_error"] < thresholds["instantaneous_oracle_error"],
        "full_chirp_velocity_matches_riccati_equation": riccati["chirp_velocity_error"] < thresholds["instantaneous_oracle_error"],
        "rotated_harmonic_oracle_generates_correlation": riccati["analytic_correlation_rate"] > 1.0e-3,
        "diagonal_width_manifold_cannot_represent_oracle_correlation": riccati["analytic_correlation_rate"] > thresholds["instantaneous_oracle_error"],
        "correlated_metric_is_positive_semidefinite": riccati["metric_minimum_eigenvalue"] > thresholds["metric_minimum_eigenvalue"],
        "correlated_metric_solve_is_accurate": riccati["metric_residual"] < thresholds["metric_residual"],
        "coarse_midpoint_matches_riccati_reference": riccati["errors"][0] < thresholds["maximum_riccati_error"],
        "riccati_error_decreases_004_to_002": riccati["errors"][1] < riccati["errors"][0],
        "riccati_error_decreases_002_to_001": riccati["errors"][2] < riccati["errors"][1],
        "riccati_convergence_is_second_order_coarse": riccati["orders"][0] > thresholds["minimum_midpoint_order"],
        "riccati_convergence_is_second_order_fine": riccati["orders"][1] > thresholds["minimum_midpoint_order"],
        "riccati_trajectory_conserves_norm_and_energy": max(riccati["maximum_norm_drift"], riccati["maximum_energy_drift"]) < thresholds["energy_drift"],
        # 61--81: arbitrary orthogonal, gauge, permutation, and reversal covariance.
        "quadratic_model_is_rotation_covariant": model_rotation_error < thresholds["covariance_error"],
        "wavefunction_is_rotation_covariant": wavefunction_rotation_error < thresholds["covariance_error"],
        "overlap_matrix_is_rotation_invariant": rotation_overlap_error < thresholds["covariance_error"],
        "hamiltonian_matrix_is_rotation_invariant": rotation_hamiltonian_error < thresholds["covariance_error"],
        "energy_is_rotation_invariant": rotation_energy_error < thresholds["covariance_error"],
        "density_is_rotation_invariant": rotation_density_error < thresholds["covariance_error"],
        "metric_is_rotation_covariant": metric_rotation_error < thresholds["covariance_error"],
        "metric_rhs_is_rotation_covariant": rhs_rotation_error < thresholds["covariance_error"],
        "tdvp_velocity_is_rotation_covariant": velocity_rotation_error < thresholds["covariance_error"],
        "implicit_midpoint_step_is_rotation_covariant": step_rotation_error < thresholds["step_covariance_error"],
        "quadratic_model_is_reflection_covariant": model_reflection_error < thresholds["covariance_error"],
        "wavefunction_is_reflection_covariant": wavefunction_reflection_error < thresholds["covariance_error"],
        "overlap_matrix_is_reflection_invariant": reflection_overlap_error < thresholds["covariance_error"],
        "hamiltonian_matrix_is_reflection_invariant": reflection_hamiltonian_error < thresholds["covariance_error"],
        "tdvp_velocity_is_reflection_covariant": velocity_reflection_error < thresholds["covariance_error"],
        "implicit_midpoint_step_is_reflection_covariant": step_reflection_error < thresholds["step_covariance_error"],
        "tdvp_velocity_is_constant_gauge_covariant": gauge_velocity_error < thresholds["covariance_error"],
        "implicit_midpoint_step_is_constant_gauge_covariant": gauge_step_error < thresholds["step_covariance_error"],
        "tdvp_velocity_is_packet_permutation_covariant": permutation_velocity_error < thresholds["covariance_error"],
        "implicit_midpoint_step_is_packet_permutation_covariant": permutation_step_error < thresholds["step_covariance_error"],
        "implicit_midpoint_is_signed_reversible": time_reversal_error < thresholds["time_reversal_error"],
        # 82--100: controlled correlated-basis lifecycle.
        "candidate_set_covers_signed_intrinsic_axes": len(candidates) == 4 * state.ndim * state.ngaussian,
        "candidate_set_covers_position_and_momentum": {item.displacement_kind for item in candidates} == {"position", "momentum"},
        "candidate_displacements_follow_principal_axes": principal_axis_error < thresholds["algebra_error"],
        "at_least_one_correlated_candidate_is_admitted": any(item.admitted for item in evaluations),
        "highest_residual_candidate_is_selected": spawn_event.selected_candidate.canonical_key() == sorted([item for item in evaluations if item.admitted], key=lambda item: (-item.residual_capture, item.candidate.canonical_key()))[0].candidate.canonical_key(),
        "spawn_adds_exactly_one_packet": spawn_event.event_kind == "spawn" and spawn_event.after.ngaussian == 2,
        "spawn_newborn_coefficient_is_exactly_zero": np.max(np.abs(spawn_event.after.coefficients[-1])) == 0.0,
        "spawn_projection_is_exact": spawn_event.projection.relative_projection_loss < thresholds["projection_loss"],
        "spawn_identity_and_age_are_stable": spawn_event.added_packet_id == "g000001" and spawn_event.packet_ages_after[-1] == 0,
        "degenerate_principal_axes_fail_closed": len(generate_correlated_spawn_candidates_v270(isotropic)) == 0,
        "near_degenerate_principal_axes_fail_closed": len(generate_correlated_spawn_candidates_v270(near_degenerate)) == 0,
        "candidate_scores_are_rotation_covariant": candidate_score_rotation_error < thresholds["candidate_score_covariance"],
        "selected_spawn_event_is_rotation_covariant": spawn_event_rotation_error < thresholds["step_covariance_error"],
        "duplicate_candidate_is_rejected_by_novelty_and_rank": (not duplicate_evaluation.admitted and duplicate_evaluation.novelty < 2.0e-12 and "rank-deficient-enlarged-basis" in duplicate_evaluation.rejection_reasons),
        "correlated_prune_passes_projection_gates": (prune_event.event_kind == "prune" and prune_event.projection.relative_projection_loss < thresholds["projection_loss"] and abs(prune_event.projection.energy_jump_hartree) < thresholds["projection_energy_jump"]),
        "correlated_merge_passes_projection_gates": (merge_event.event_kind == "merge" and merge_event.projection.relative_projection_loss < thresholds["projection_loss"] and abs(merge_event.projection.energy_jump_hartree) < thresholds["projection_energy_jump"]),
        "newborn_full_shape_block_starts_dormant": dormant_mask.tolist() == [True, False],
        "newborn_full_shape_activation_requires_metric_gate": activated_mask.tolist() == [True, True],
        "dormant_full_matrix_shape_is_exactly_frozen": dormant_shape_drift < thresholds["dormant_shape_drift"],
    }
    if len(checks) != 100:
        raise AssertionError(f"v0.27.0 must define exactly 100 validation gates, found {len(checks)}.")
    evidence = CorrelatedValidationEvidenceV270(
        metrics=metrics,
        thresholds=thresholds,
        checks={name: bool(value) for name, value in checks.items()},
        trajectory_fingerprint=riccati["trajectory_fingerprint"],
        lifecycle_fingerprint=_sha256_v270(spawn_event.as_dict()),
    )
    return evidence.validate()


def save_correlated_validation_evidence_v270(path):
    evidence = run_correlated_validation_evidence_v270()
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(evidence.as_dict(), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return evidence
