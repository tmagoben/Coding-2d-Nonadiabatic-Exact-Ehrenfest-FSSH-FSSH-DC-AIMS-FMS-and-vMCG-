"""Adaptive-width one-dimensional multi-Gaussian TDVP for v0.25.2.

Each normalized packet carries a positive real width and a real quadratic chirp,

    g_I(x) = (alpha_I/pi)**(1/4)
             exp[-alpha_I (x-q_I)**2/2
                 + i beta_I (x-q_I)**2/2
                 + i p_I (x-q_I)].

The propagated real coordinates are coefficient real/imaginary parts, centers,
momenta, logarithmic widths, and chirps.  Logarithmic widths make positivity
structural.  Exact complex Gaussian moments through degree four build the
McLachlan metric and RHS.  The metric uses the inherited compatible-null full-SVD
solve and the complete vector field is advanced by fully implicit midpoint.

The released scope remains one nuclear coordinate, a fixed complete electronic
frame, constant positive mass, and a Hermitian matrix potential through degree two.
Spawning, pruning, multidimensional width matrices, coordinate-dependent electronic
frames, and real molecular-SOC trajectories remain closed.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np
from scipy.optimize import root

from .multigaussian_tdvp_v251 import (
    MetricSolveReceiptV251,
    QuadraticSpinHamiltonianV251,
    quadratic_spin_hamiltonian_from_provider_v251,
    solve_variational_metric_v251,
)


ADAPTIVE_MULTIGAUSSIAN_TDVP_SCHEMA_V252 = (
    "gnd-adaptive-width-multigaussian-tdvp-trajectory-v0.25.2"
)
ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252 = (
    "one-dimensional coupled normalized Gaussians with logarithmic widths, "
    "quadratic chirps, and complete fixed-frame spinors"
)
VARIATIONAL_PRINCIPLE_V252 = "real-parameter McLachlan time-dependent variation"
VARIATIONAL_INTEGRATOR_V252 = (
    "fully implicit midpoint applied to the adaptive-width SVD-pseudoinverse "
    "TDVP vector field"
)
VARIATIONAL_METRIC_SOLVER_V252 = (
    "full SVD with relative/absolute rank cutoff and compatible-null-space audit"
)
WIDTH_COORDINATES_V252 = (
    "eta=log(alpha) for positive width alpha, paired with real quadratic chirp beta"
)
POTENTIAL_CONTRACT_V252 = (
    "fixed-frame one-dimensional Hermitian quadratic complete-spinor Hamiltonian"
)

# The v0.25.2 adaptive layer deliberately inherits the already audited quadratic
# Hamiltonian object while exposing a release-local public name.
QuadraticSpinHamiltonianV252 = QuadraticSpinHamiltonianV251


def _canonical_v252(value):
    if isinstance(value, np.generic):
        return _canonical_v252(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_v252(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_v252(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v252(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("v0.25.2 canonical data cannot contain non-finite values.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.25.2 canonical value: {type(value).__name__}")


def _sha256_v252(value):
    payload = json.dumps(
        _canonical_v252(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_v252(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_norm_v252(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def _hermiticity_residual_v252(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    return _scaled_norm_v252(matrix, matrix.conj().T)


def _unitary_v252(matrix, *, tolerance):
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("electronic gauge matrix must be nonempty and square.")
    identity = np.eye(matrix.shape[0], dtype=complex)
    if _scaled_norm_v252(matrix.conj().T @ matrix, identity) > float(tolerance):
        raise ValueError("electronic gauge matrix must be unitary.")
    return matrix


@dataclass(frozen=True)
class AdaptiveVariationalSettingsV252:
    """Frozen algorithms and width-domain gates for the v0.25.2 solver."""

    metric_relative_cutoff: float = 1.0e-10
    metric_absolute_cutoff: float = 1.0e-12
    maximum_retained_condition_number: float = 1.0e10
    null_rhs_relative_tolerance: float = 2.0e-9
    linear_residual_relative_tolerance: float = 2.0e-9
    nonlinear_residual_tolerance: float = 2.0e-10
    nonlinear_xtol: float = 1.0e-10
    nonlinear_max_function_evaluations: int = 800
    structural_tolerance: float = 2.0e-10
    maximum_step_norm_drift: float = 3.0e-8
    minimum_width: float = 1.0e-8
    maximum_width: float = 1.0e8
    maximum_absolute_chirp: float = 1.0e8
    maximum_step_log_width_change: float = 0.5
    variational_principle: str = VARIATIONAL_PRINCIPLE_V252
    integrator: str = VARIATIONAL_INTEGRATOR_V252
    metric_solver: str = VARIATIONAL_METRIC_SOLVER_V252
    nonlinear_solver: str = "scipy.optimize.root-hybr"
    width_coordinates: str = WIDTH_COORDINATES_V252
    allow_compatible_rank_deficiency: bool = True
    adaptive_gaussian_widths: bool = True
    spawning: bool = False
    pruning: bool = False
    coordinate_dependent_electronic_frame: bool = False
    multidimensional_nuclear_motion: bool = False
    full_width_matrices: bool = False
    real_molecular_soc_provider: bool = False

    def validate(self):
        for name in (
            "metric_relative_cutoff",
            "metric_absolute_cutoff",
            "maximum_retained_condition_number",
            "null_rhs_relative_tolerance",
            "linear_residual_relative_tolerance",
            "nonlinear_residual_tolerance",
            "nonlinear_xtol",
            "structural_tolerance",
            "maximum_step_norm_drift",
            "minimum_width",
            "maximum_width",
            "maximum_absolute_chirp",
            "maximum_step_log_width_change",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if float(self.metric_relative_cutoff) >= 1.0:
            raise ValueError("metric_relative_cutoff must be smaller than one.")
        if float(self.maximum_retained_condition_number) < 1.0:
            raise ValueError("maximum_retained_condition_number must be at least one.")
        if float(self.maximum_width) <= float(self.minimum_width):
            raise ValueError("maximum_width must exceed minimum_width.")
        if (
            isinstance(self.nonlinear_max_function_evaluations, (bool, np.bool_))
            or not isinstance(
                self.nonlinear_max_function_evaluations, (int, np.integer)
            )
            or int(self.nonlinear_max_function_evaluations) < 1
        ):
            raise ValueError(
                "nonlinear_max_function_evaluations must be a positive integer."
            )
        for name in (
            "allow_compatible_rank_deficiency",
            "adaptive_gaussian_widths",
            "spawning",
            "pruning",
            "coordinate_dependent_electronic_frame",
            "multidimensional_nuclear_motion",
            "full_width_matrices",
            "real_molecular_soc_provider",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if not self.allow_compatible_rank_deficiency:
            raise ValueError(
                "v0.25.2 freezes compatible SVD null-space handling as enabled."
            )
        if not self.adaptive_gaussian_widths:
            raise ValueError("v0.25.2 freezes adaptive Gaussian widths as enabled.")
        closed = {
            "spawning": self.spawning,
            "pruning": self.pruning,
            "coordinate_dependent_electronic_frame": (
                self.coordinate_dependent_electronic_frame
            ),
            "multidimensional_nuclear_motion": self.multidimensional_nuclear_motion,
            "full_width_matrices": self.full_width_matrices,
            "real_molecular_soc_provider": self.real_molecular_soc_provider,
        }
        requested = [name for name, enabled in closed.items() if enabled]
        if requested:
            raise ValueError(
                "v0.25.2 does not admit: " + ", ".join(sorted(requested)) + "."
            )
        if self.variational_principle != VARIATIONAL_PRINCIPLE_V252:
            raise ValueError("the v0.25.2 variational principle is frozen.")
        if self.integrator != VARIATIONAL_INTEGRATOR_V252:
            raise ValueError("the v0.25.2 implicit integrator is frozen.")
        if self.metric_solver != VARIATIONAL_METRIC_SOLVER_V252:
            raise ValueError("the v0.25.2 SVD metric solver is frozen.")
        if self.nonlinear_solver != "scipy.optimize.root-hybr":
            raise ValueError("the v0.25.2 nonlinear solver is frozen.")
        if self.width_coordinates != WIDTH_COORDINATES_V252:
            raise ValueError("the v0.25.2 width coordinates are frozen.")
        return self

    def as_dict(self):
        self.validate()
        return _canonical_v252(asdict(self))


def quadratic_spin_hamiltonian_from_provider_v252(provider, *, tolerance=3.0e-11):
    """Reuse the independently audited v0.25.1 quadratic provider intake."""

    return quadratic_spin_hamiltonian_from_provider_v251(
        provider, tolerance=tolerance
    )


def _cross_moments_v252(
    qi,
    pi,
    alphai,
    betai,
    qj,
    pj,
    alphaj,
    betaj,
    maximum_order=4,
):
    """Return exact <g_i|x^n|g_j> moments through degree four."""

    if (
        isinstance(maximum_order, (bool, np.bool_))
        or not isinstance(maximum_order, (int, np.integer))
        or not 0 <= int(maximum_order) <= 4
    ):
        raise ValueError("v0.25.2 cross moments are implemented through order four.")
    values = np.asarray(
        [qi, pi, alphai, betai, qj, pj, alphaj, betaj], dtype=float
    )
    if not np.all(np.isfinite(values)):
        raise ValueError("adaptive Gaussian parameters must be finite.")
    qi, pi, alphai, betai, qj, pj, alphaj, betaj = values.tolist()
    if alphai <= 0.0 or alphaj <= 0.0:
        raise ValueError("adaptive Gaussian widths must be positive.")

    zi_bra = alphai + 1.0j * betai
    zj_ket = alphaj - 1.0j * betaj
    combined = zi_bra + zj_ket
    linear = zi_bra * qi + zj_ket * qj + 1.0j * (pj - pi)
    constant = (
        -0.5 * zi_bra * qi**2
        - 0.5 * zj_ket * qj**2
        + 1.0j * pi * qi
        - 1.0j * pj * qj
    )
    log_prefactor = (
        0.25 * np.log(alphai)
        + 0.25 * np.log(alphaj)
        + 0.5 * np.log(2.0)
        - 0.5 * np.log(combined)
    )
    overlap = np.exp(
        log_prefactor + constant + 0.5 * linear * linear / combined
    )
    mean = linear / combined
    variance = 1.0 / combined
    moments = [overlap]
    if maximum_order >= 1:
        moments.append(overlap * mean)
    if maximum_order >= 2:
        moments.append(overlap * (mean**2 + variance))
    if maximum_order >= 3:
        moments.append(overlap * (mean**3 + 3.0 * mean * variance))
    if maximum_order >= 4:
        moments.append(
            overlap * (mean**4 + 6.0 * mean**2 * variance + 3.0 * variance**2)
        )
    result = np.asarray(moments, dtype=complex)
    if not np.all(np.isfinite(result)):
        raise ValueError("adaptive Gaussian cross moments are non-finite.")
    return result


def _integrate_polynomial_v252(coefficients, moments):
    coefficients = np.asarray(coefficients, dtype=complex)
    moments = np.asarray(moments, dtype=complex)
    if coefficients.ndim != 1 or len(coefficients) > len(moments):
        raise ValueError("polynomial degree exceeds the available Gaussian moments.")
    return complex(np.dot(coefficients, moments[: len(coefficients)]))


def _shifted_quadratic_to_x_v252(constant, linear, quadratic, center):
    center = float(center)
    return np.asarray(
        [
            constant - linear * center + quadratic * center**2,
            linear - 2.0 * quadratic * center,
            quadratic,
        ],
        dtype=complex,
    )


def _kinetic_polynomial_v252(q, p, width, chirp, mass):
    q = float(q)
    p = float(p)
    width = float(width)
    chirp = float(chirp)
    mass = float(mass)
    z = width - 1.0j * chirp
    return _shifted_quadratic_to_x_v252(
        (p**2 + z) / (2.0 * mass),
        1.0j * z * p / mass,
        -(z**2) / (2.0 * mass),
        q,
    )


@dataclass(frozen=True)
class ThawedGaussianSpinorStateV252:
    q: np.ndarray
    p: np.ndarray
    widths: np.ndarray
    chirps: np.ndarray
    coefficients: np.ndarray
    time_au: float = 0.0

    def __post_init__(self):
        for name in ("q", "p", "widths", "chirps"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=float).copy()
            )
        object.__setattr__(
            self,
            "coefficients",
            np.asarray(self.coefficients, dtype=complex).copy(),
        )

    @property
    def ngaussian(self):
        return int(len(self.q)) if self.q.ndim == 1 else 0

    @property
    def nstate(self):
        return int(self.coefficients.shape[1]) if self.coefficients.ndim == 2 else 0

    @property
    def parameter_count(self):
        return 2 * self.ngaussian * self.nstate + 4 * self.ngaussian

    @property
    def log_widths(self):
        if self.widths.ndim != 1 or np.any(self.widths <= 0.0):
            raise ValueError("cannot form logarithmic coordinates from invalid widths.")
        return np.log(self.widths)

    def nuclear_overlap_matrix(self):
        if (
            self.q.ndim != 1
            or len(self.q) < 1
            or self.p.shape != self.q.shape
            or self.widths.shape != self.q.shape
            or self.chirps.shape != self.q.shape
            or np.any(self.widths <= 0.0)
        ):
            raise ValueError("cannot build overlap from invalid adaptive packets.")
        overlap = np.zeros((self.ngaussian, self.ngaussian), dtype=complex)
        for i in range(self.ngaussian):
            for j in range(self.ngaussian):
                overlap[i, j] = _cross_moments_v252(
                    self.q[i],
                    self.p[i],
                    self.widths[i],
                    self.chirps[i],
                    self.q[j],
                    self.p[j],
                    self.widths[j],
                    self.chirps[j],
                    maximum_order=0,
                )[0]
        return overlap

    @property
    def generalized_norm(self):
        if self.coefficients.ndim != 2 or self.q.ndim != 1:
            return float("nan")
        overlap = self.nuclear_overlap_matrix()
        return float(
            np.real(
                np.einsum(
                    "ia,ij,ja->",
                    np.conj(self.coefficients),
                    overlap,
                    self.coefficients,
                )
            )
        )

    def validate(self, *, require_normalized=True, tolerance=2.0e-10):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("state tolerance must be finite and positive.")
        if self.q.ndim != 1 or len(self.q) < 1:
            raise ValueError("adaptive Gaussian coordinates must be a nonempty vector.")
        if any(
            value.shape != self.q.shape
            for value in (self.p, self.widths, self.chirps)
        ):
            raise ValueError("adaptive Gaussian q, p, width, and chirp shapes differ.")
        if self.coefficients.ndim != 2 or self.coefficients.shape[0] != len(self.q):
            raise ValueError("coefficient array must have shape (n_gaussian,n_state).")
        if self.coefficients.shape[1] < 1:
            raise ValueError("adaptive Gaussian state requires an electronic state.")
        if not all(
            np.all(np.isfinite(value))
            for value in (
                self.q,
                self.p,
                self.widths,
                self.chirps,
                self.coefficients,
            )
        ):
            raise ValueError("adaptive Gaussian state contains non-finite values.")
        if np.any(self.widths <= 0.0):
            raise ValueError("all adaptive Gaussian widths must be positive.")
        if not np.isfinite(float(self.time_au)):
            raise ValueError("adaptive Gaussian time must be finite.")
        norm = self.generalized_norm
        if not np.isfinite(norm) or norm <= tolerance:
            raise ValueError("adaptive Gaussian wavefunction norm must be positive.")
        if require_normalized and abs(norm - 1.0) > tolerance:
            raise ValueError(
                "adaptive Gaussian wavefunction must have unit generalized norm."
            )
        return self

    def normalized(self):
        self.validate(require_normalized=False)
        return ThawedGaussianSpinorStateV252(
            q=self.q,
            p=self.p,
            widths=self.widths,
            chirps=self.chirps,
            coefficients=self.coefficients / np.sqrt(self.generalized_norm),
            time_au=self.time_au,
        ).validate()

    def permuted(self, order):
        order = np.asarray(order)
        if (
            order.shape != (self.ngaussian,)
            or not np.issubdtype(order.dtype, np.integer)
            or set(order.tolist()) != set(range(self.ngaussian))
        ):
            raise ValueError("Gaussian permutation must contain every packet once.")
        return ThawedGaussianSpinorStateV252(
            q=self.q[order],
            p=self.p[order],
            widths=self.widths[order],
            chirps=self.chirps[order],
            coefficients=self.coefficients[order],
            time_au=self.time_au,
        ).validate(tolerance=2.0e-9)

    def gauge_transformed(self, unitary):
        unitary = _unitary_v252(unitary, tolerance=2.0e-11)
        if unitary.shape != (self.nstate, self.nstate):
            raise ValueError("electronic gauge dimension disagrees with the state.")
        return ThawedGaussianSpinorStateV252(
            q=self.q,
            p=self.p,
            widths=self.widths,
            chirps=self.chirps,
            coefficients=self.coefficients @ unitary.conj(),
            time_au=self.time_au,
        ).validate(tolerance=2.0e-9)

    def as_dict(self):
        self.validate(require_normalized=False)
        return {
            "q": self.q.tolist(),
            "p": self.p.tolist(),
            "widths": self.widths.tolist(),
            "log_widths": self.log_widths.tolist(),
            "chirps": self.chirps.tolist(),
            "coefficients": _complex_pairs_v252(self.coefficients),
            "time_au": float(self.time_au),
            "generalized_norm": self.generalized_norm,
            "n_gaussian": self.ngaussian,
            "n_state": self.nstate,
            "parameter_count": self.parameter_count,
        }


def _validate_width_domain_v252(state, settings):
    state = state.validate(require_normalized=False)
    settings = settings.validate()
    if np.min(state.widths) < settings.minimum_width:
        raise ValueError("adaptive Gaussian width fell below the configured minimum.")
    if np.max(state.widths) > settings.maximum_width:
        raise ValueError("adaptive Gaussian width exceeded the configured maximum.")
    if np.max(np.abs(state.chirps)) > settings.maximum_absolute_chirp:
        raise ValueError("adaptive Gaussian chirp exceeded the configured maximum.")
    return state


def pack_adaptive_variational_parameters_v252(state):
    state = state.validate(require_normalized=False)
    return np.concatenate(
        (
            state.coefficients.real.reshape(-1),
            state.coefficients.imag.reshape(-1),
            state.q,
            state.p,
            state.log_widths,
            state.chirps,
        )
    )


def state_from_adaptive_variational_parameters_v252(
    parameters,
    *,
    ngaussian,
    nstate,
    time_au,
):
    parameters = np.asarray(parameters, dtype=float)
    if (
        isinstance(ngaussian, (bool, np.bool_))
        or not isinstance(ngaussian, (int, np.integer))
        or int(ngaussian) < 1
    ):
        raise ValueError("ngaussian must be a positive integer.")
    if (
        isinstance(nstate, (bool, np.bool_))
        or not isinstance(nstate, (int, np.integer))
        or int(nstate) < 1
    ):
        raise ValueError("nstate must be a positive integer.")
    ngaussian = int(ngaussian)
    nstate = int(nstate)
    expected = 2 * ngaussian * nstate + 4 * ngaussian
    if parameters.shape != (expected,) or not np.all(np.isfinite(parameters)):
        raise ValueError("adaptive variational parameter vector is invalid.")
    block = ngaussian * nstate
    coefficients = (
        parameters[:block] + 1.0j * parameters[block : 2 * block]
    ).reshape(ngaussian, nstate)
    cursor = 2 * block
    q = parameters[cursor : cursor + ngaussian]
    cursor += ngaussian
    p = parameters[cursor : cursor + ngaussian]
    cursor += ngaussian
    log_widths = parameters[cursor : cursor + ngaussian]
    cursor += ngaussian
    if np.max(np.abs(log_widths)) > 700.0:
        raise ValueError("logarithmic Gaussian width is outside floating-point range.")
    widths = np.exp(log_widths)
    chirps = parameters[cursor : cursor + ngaussian]
    return ThawedGaussianSpinorStateV252(
        q=q,
        p=p,
        widths=widths,
        chirps=chirps,
        coefficients=coefficients,
        time_au=float(time_au),
    ).validate(require_normalized=False)


def build_adaptive_gaussian_spinor_matrices_v252(state, model):
    """Return exact overlap/H matrices for chirped adaptive packets."""

    state = state.validate(require_normalized=False)
    model = model.validate()
    if state.nstate != model.nstate:
        raise ValueError("state and Hamiltonian electronic dimensions differ.")
    dimension = state.ngaussian * state.nstate
    overlap = np.zeros((dimension, dimension), dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(state.nstate, dtype=complex)
    for i in range(state.ngaussian):
        si = slice(i * state.nstate, (i + 1) * state.nstate)
        for j in range(state.ngaussian):
            sj = slice(j * state.nstate, (j + 1) * state.nstate)
            moments = _cross_moments_v252(
                state.q[i],
                state.p[i],
                state.widths[i],
                state.chirps[i],
                state.q[j],
                state.p[j],
                state.widths[j],
                state.chirps[j],
                maximum_order=2,
            )
            kinetic = _kinetic_polynomial_v252(
                state.q[j],
                state.p[j],
                state.widths[j],
                state.chirps[j],
                model.mass_au,
            )
            overlap[si, sj] = moments[0] * identity
            hamiltonian[si, sj] = (
                _integrate_polynomial_v252(kinetic, moments) * identity
                + moments[0] * model.H0
                + moments[1] * model.H1
                + moments[2] * model.H2
            )
    if _hermiticity_residual_v252(overlap) > 2.0e-10:
        raise ValueError("adaptive Gaussian overlap matrix is not Hermitian.")
    if _hermiticity_residual_v252(hamiltonian) > 2.0e-10:
        raise ValueError("adaptive Gaussian Hamiltonian matrix is not Hermitian.")
    return overlap, hamiltonian


def adaptive_variational_energy_v252(state, model):
    state = state.validate(require_normalized=False)
    overlap, hamiltonian = build_adaptive_gaussian_spinor_matrices_v252(state, model)
    coefficients = state.coefficients.reshape(-1)
    norm = float(np.real(np.vdot(coefficients, overlap @ coefficients)))
    if norm <= 0.0:
        raise ValueError("cannot evaluate energy of a zero-norm adaptive state.")
    value = np.vdot(coefficients, hamiltonian @ coefficients) / norm
    if abs(float(np.imag(value))) > 2.0e-10:
        raise ValueError("adaptive variational energy has an imaginary component.")
    return float(np.real(value))


def _tangent_terms_v252(state):
    """Return (packet, electronic vector, x-polynomial) in pack order."""

    terms = []
    nstate = state.nstate
    for i in range(state.ngaussian):
        for electronic in range(nstate):
            vector = np.zeros(nstate, dtype=complex)
            vector[electronic] = 1.0
            terms.append((i, vector, np.asarray([1.0 + 0.0j])))
    for i in range(state.ngaussian):
        for electronic in range(nstate):
            vector = np.zeros(nstate, dtype=complex)
            vector[electronic] = 1.0j
            terms.append((i, vector, np.asarray([1.0 + 0.0j])))
    for i in range(state.ngaussian):
        z = state.widths[i] - 1.0j * state.chirps[i]
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray([-z * state.q[i] - 1.0j * state.p[i], z]),
            )
        )
    for i in range(state.ngaussian):
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray([-1.0j * state.q[i], 1.0j]),
            )
        )
    for i in range(state.ngaussian):
        q = state.q[i]
        width = state.widths[i]
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray(
                    [0.25 - 0.5 * width * q**2, width * q, -0.5 * width],
                    dtype=complex,
                ),
            )
        )
    for i in range(state.ngaussian):
        q = state.q[i]
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray([0.5j * q**2, -1.0j * q, 0.5j]),
            )
        )
    if len(terms) != state.parameter_count:
        raise AssertionError("v0.25.2 tangent count disagrees with parameter layout.")
    return terms


@dataclass(frozen=True)
class AdaptiveVariationalMetricSystemV252:
    metric: np.ndarray
    rhs: np.ndarray
    velocity: np.ndarray
    solve_receipt: MetricSolveReceiptV251
    generalized_norm: float
    energy_hartree: float
    settings: AdaptiveVariationalSettingsV252

    def __post_init__(self):
        for name in ("metric", "rhs", "velocity"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=float).copy()
            )

    def validate(self):
        settings = self.settings.validate()
        if self.metric.ndim != 2 or self.metric.shape[0] != self.metric.shape[1]:
            raise ValueError("stored adaptive variational metric must be square.")
        count = self.metric.shape[0]
        if self.rhs.shape != (count,) or self.velocity.shape != (count,):
            raise ValueError("stored adaptive metric vectors have incompatible shape.")
        if not all(
            np.all(np.isfinite(value))
            for value in (self.metric, self.rhs, self.velocity)
        ):
            raise ValueError("stored adaptive variational metric is non-finite.")
        expected_velocity, expected_receipt = solve_variational_metric_v251(
            self.metric, self.rhs, settings=settings
        )
        if _scaled_norm_v252(self.velocity, expected_velocity) > settings.structural_tolerance:
            raise ValueError("stored adaptive TDVP velocity disagrees with SVD solve.")
        if _scaled_norm_v252(
            self.solve_receipt.singular_values,
            expected_receipt.singular_values,
        ) > settings.structural_tolerance:
            raise ValueError("stored adaptive metric spectrum is inconsistent.")
        for name in (
            "cutoff",
            "retained_condition_number",
            "rhs_norm",
            "null_rhs_norm",
            "null_rhs_relative",
            "linear_residual_norm",
            "linear_residual_relative",
            "velocity_norm",
        ):
            if abs(
                float(getattr(self.solve_receipt, name))
                - float(getattr(expected_receipt, name))
            ) > settings.structural_tolerance:
                raise ValueError(f"stored adaptive metric receipt {name} is inconsistent.")
        if (
            self.solve_receipt.rank != expected_receipt.rank
            or self.solve_receipt.nullity != expected_receipt.nullity
        ):
            raise ValueError("stored adaptive metric rank/nullity are inconsistent.")
        if not np.isfinite(float(self.generalized_norm)) or float(
            self.generalized_norm
        ) <= 0.0:
            raise ValueError("stored adaptive generalized norm must be positive.")
        if not np.isfinite(float(self.energy_hartree)):
            raise ValueError("stored adaptive energy must be finite.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "metric": self.metric.tolist(),
            "rhs": self.rhs.tolist(),
            "velocity": self.velocity.tolist(),
            "solve_receipt": self.solve_receipt.as_dict(),
            "generalized_norm": float(self.generalized_norm),
            "energy_hartree": float(self.energy_hartree),
            "settings": self.settings.as_dict(),
        }


def build_adaptive_variational_metric_system_v252(
    state,
    model,
    *,
    settings=AdaptiveVariationalSettingsV252(),
):
    """Build and solve the exact adaptive-width McLachlan metric system."""

    settings = settings.validate()
    state = _validate_width_domain_v252(state, settings)
    model = model.validate()
    if state.nstate != model.nstate:
        raise ValueError("adaptive state and model electronic dimensions differ.")
    terms = _tangent_terms_v252(state)
    count = len(terms)
    metric = np.zeros((count, count), dtype=float)
    rhs = np.zeros(count, dtype=float)
    moment_cache = {}
    for i in range(state.ngaussian):
        for j in range(state.ngaussian):
            moment_cache[(i, j)] = _cross_moments_v252(
                state.q[i],
                state.p[i],
                state.widths[i],
                state.chirps[i],
                state.q[j],
                state.p[j],
                state.widths[j],
                state.chirps[j],
                maximum_order=4,
            )

    for mu, (i, vector_i, polynomial_i) in enumerate(terms):
        conjugate_i = np.conj(polynomial_i)
        for nu, (j, vector_j, polynomial_j) in enumerate(terms):
            polynomial = np.convolve(conjugate_i, polynomial_j)
            nuclear = _integrate_polynomial_v252(
                polynomial, moment_cache[(i, j)]
            )
            metric[mu, nu] = float(np.real(np.vdot(vector_i, vector_j) * nuclear))

        projection = 0.0 + 0.0j
        for j in range(state.ngaussian):
            kinetic = _kinetic_polynomial_v252(
                state.q[j],
                state.p[j],
                state.widths[j],
                state.chirps[j],
                model.mass_au,
            )
            matrix_polynomial = (
                kinetic[0] * np.eye(model.nstate, dtype=complex) + model.H0,
                kinetic[1] * np.eye(model.nstate, dtype=complex) + model.H1,
                kinetic[2] * np.eye(model.nstate, dtype=complex) + model.H2,
            )
            moments = moment_cache[(i, j)]
            for bra_degree, bra_value in enumerate(conjugate_i):
                for operator_degree, operator_matrix in enumerate(matrix_polynomial):
                    projection += (
                        bra_value
                        * moments[bra_degree + operator_degree]
                        * np.vdot(
                            vector_i,
                            operator_matrix @ state.coefficients[j],
                        )
                    )
        rhs[mu] = float(np.imag(projection))

    metric = 0.5 * (metric + metric.T)
    velocity, solve_receipt = solve_variational_metric_v251(
        metric, rhs, settings=settings
    )
    return AdaptiveVariationalMetricSystemV252(
        metric=metric,
        rhs=rhs,
        velocity=velocity,
        solve_receipt=solve_receipt,
        generalized_norm=state.generalized_norm,
        energy_hartree=adaptive_variational_energy_v252(state, model),
        settings=settings,
    ).validate()


@dataclass(frozen=True)
class AdaptiveImplicitMidpointTDVPStepV252:
    start: ThawedGaussianSpinorStateV252
    end: ThawedGaussianSpinorStateV252
    dt_au: float
    model: QuadraticSpinHamiltonianV251
    settings: AdaptiveVariationalSettingsV252
    midpoint_parameters: np.ndarray
    midpoint_system: AdaptiveVariationalMetricSystemV252
    nonlinear_success: bool
    nonlinear_status: int
    nonlinear_message: str
    nonlinear_function_evaluations: int
    nonlinear_jacobian_evaluations: int | None
    nonlinear_residual: np.ndarray
    nonlinear_residual_norm: float
    predictor_residual_norm: float
    start_norm: float
    end_norm: float
    start_energy_hartree: float
    end_energy_hartree: float

    def __post_init__(self):
        for name in ("midpoint_parameters", "nonlinear_residual"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=float).copy()
            )

    @property
    def norm_change(self):
        return float(self.end_norm - self.start_norm)

    @property
    def energy_change_hartree(self):
        return float(self.end_energy_hartree - self.start_energy_hartree)

    @property
    def maximum_log_width_change(self):
        return float(np.max(np.abs(self.end.log_widths - self.start.log_widths)))

    def validate(self):
        settings = self.settings.validate()
        model = self.model.validate()
        state_tolerance = max(
            settings.structural_tolerance, settings.maximum_step_norm_drift
        )
        self.start.validate(require_normalized=True, tolerance=state_tolerance)
        self.end.validate(require_normalized=False, tolerance=state_tolerance)
        _validate_width_domain_v252(self.start, settings)
        _validate_width_domain_v252(self.end, settings)
        if self.start.nstate != model.nstate or self.end.nstate != model.nstate:
            raise ValueError("adaptive step state/model dimensions differ.")
        if self.start.ngaussian != self.end.ngaussian:
            raise ValueError("Gaussian count changed inside a v0.25.2 step.")
        self.end.validate(require_normalized=True, tolerance=state_tolerance)
        if self.maximum_log_width_change > settings.maximum_step_log_width_change:
            raise ValueError("adaptive width changed too much inside one step.")
        dt = float(self.dt_au)
        if not np.isfinite(dt) or dt == 0.0:
            raise ValueError("adaptive implicit step dt must be finite and nonzero.")
        if abs(self.end.time_au - (self.start.time_au + dt)) > settings.structural_tolerance:
            raise ValueError("adaptive implicit endpoint time disagrees with dt.")
        theta_start = pack_adaptive_variational_parameters_v252(self.start)
        theta_end = pack_adaptive_variational_parameters_v252(self.end)
        expected_midpoint = 0.5 * (theta_start + theta_end)
        if _scaled_norm_v252(
            self.midpoint_parameters, expected_midpoint
        ) > settings.structural_tolerance:
            raise ValueError("stored adaptive midpoint parameters are inconsistent.")
        midpoint_state = state_from_adaptive_variational_parameters_v252(
            expected_midpoint,
            ngaussian=self.start.ngaussian,
            nstate=self.start.nstate,
            time_au=self.start.time_au + 0.5 * dt,
        )
        expected_system = build_adaptive_variational_metric_system_v252(
            midpoint_state, model, settings=settings
        )
        self.midpoint_system.validate()
        if self.midpoint_system.settings.as_dict() != settings.as_dict():
            raise ValueError("stored adaptive midpoint settings are inconsistent.")
        for name in ("metric", "rhs", "velocity"):
            if _scaled_norm_v252(
                getattr(self.midpoint_system, name), getattr(expected_system, name)
            ) > settings.structural_tolerance:
                raise ValueError(f"stored adaptive midpoint {name} is inconsistent.")
        if _scaled_norm_v252(
            self.midpoint_system.solve_receipt.singular_values,
            expected_system.solve_receipt.singular_values,
        ) > settings.structural_tolerance:
            raise ValueError("stored adaptive midpoint spectrum is inconsistent.")
        expected_residual = theta_end - theta_start - dt * expected_system.velocity
        if self.nonlinear_residual.shape != expected_residual.shape or _scaled_norm_v252(
            self.nonlinear_residual, expected_residual
        ) > settings.structural_tolerance:
            raise ValueError("stored adaptive nonlinear residual is inconsistent.")
        expected_residual_norm = float(np.linalg.norm(expected_residual))
        if abs(
            float(self.nonlinear_residual_norm) - expected_residual_norm
        ) > settings.structural_tolerance:
            raise ValueError("stored adaptive nonlinear residual norm is inconsistent.")
        if expected_residual_norm > settings.nonlinear_residual_tolerance:
            raise ValueError("adaptive midpoint nonlinear residual exceeds tolerance.")
        if type(self.nonlinear_success) is not bool or not self.nonlinear_success:
            raise ValueError("adaptive midpoint nonlinear solver did not report success.")
        if (
            isinstance(self.nonlinear_status, (bool, np.bool_))
            or not isinstance(self.nonlinear_status, (int, np.integer))
            or int(self.nonlinear_status) <= 0
        ):
            raise ValueError("adaptive midpoint nonlinear status is not successful.")
        if not isinstance(self.nonlinear_message, str) or not self.nonlinear_message:
            raise ValueError("adaptive midpoint nonlinear message is missing.")
        if (
            isinstance(self.nonlinear_function_evaluations, (bool, np.bool_))
            or not isinstance(
                self.nonlinear_function_evaluations, (int, np.integer)
            )
            or int(self.nonlinear_function_evaluations) < 1
        ):
            raise ValueError("adaptive midpoint function-evaluation count is invalid.")
        if self.nonlinear_jacobian_evaluations is not None:
            if (
                isinstance(self.nonlinear_jacobian_evaluations, (bool, np.bool_))
                or not isinstance(
                    self.nonlinear_jacobian_evaluations, (int, np.integer)
                )
                or int(self.nonlinear_jacobian_evaluations) < 0
            ):
                raise ValueError(
                    "adaptive midpoint Jacobian-evaluation count is invalid."
                )
        if not np.isfinite(float(self.predictor_residual_norm)) or float(
            self.predictor_residual_norm
        ) < 0.0:
            raise ValueError("adaptive midpoint predictor residual is invalid.")
        expected_values = {
            "start_norm": self.start.generalized_norm,
            "end_norm": self.end.generalized_norm,
            "start_energy_hartree": adaptive_variational_energy_v252(
                self.start, model
            ),
            "end_energy_hartree": adaptive_variational_energy_v252(self.end, model),
        }
        for name, expected in expected_values.items():
            if abs(float(getattr(self, name)) - float(expected)) > settings.structural_tolerance:
                raise ValueError(f"stored adaptive {name} is inconsistent.")
        if abs(self.norm_change) > settings.maximum_step_norm_drift:
            raise ValueError("adaptive TDVP step norm drift exceeds the release gate.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "start": self.start.as_dict(),
            "end": self.end.as_dict(),
            "dt_au": float(self.dt_au),
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "midpoint_parameters": self.midpoint_parameters.tolist(),
            "midpoint_system": self.midpoint_system.as_dict(),
            "nonlinear": {
                "success": self.nonlinear_success,
                "status": int(self.nonlinear_status),
                "message": self.nonlinear_message,
                "function_evaluations": int(self.nonlinear_function_evaluations),
                "jacobian_evaluations": (
                    None
                    if self.nonlinear_jacobian_evaluations is None
                    else int(self.nonlinear_jacobian_evaluations)
                ),
                "residual": self.nonlinear_residual.tolist(),
                "residual_norm": float(self.nonlinear_residual_norm),
                "predictor_residual_norm": float(self.predictor_residual_norm),
            },
            "start_norm": float(self.start_norm),
            "end_norm": float(self.end_norm),
            "norm_change": self.norm_change,
            "start_energy_hartree": float(self.start_energy_hartree),
            "end_energy_hartree": float(self.end_energy_hartree),
            "energy_change_hartree": self.energy_change_hartree,
            "maximum_log_width_change": self.maximum_log_width_change,
        }


def adaptive_implicit_midpoint_tdvp_step_v252(
    state,
    model,
    dt_au,
    *,
    settings=AdaptiveVariationalSettingsV252(),
):
    """Advance one signed, fully implicit adaptive-width TDVP step."""

    settings = settings.validate()
    model = model.validate()
    state = state.validate(
        require_normalized=True,
        tolerance=max(
            settings.structural_tolerance, settings.maximum_step_norm_drift
        ),
    )
    _validate_width_domain_v252(state, settings)
    if state.nstate != model.nstate:
        raise ValueError("adaptive state and model electronic dimensions differ.")
    dt = float(dt_au)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError("adaptive implicit step dt must be finite and nonzero.")
    theta_start = pack_adaptive_variational_parameters_v252(state)
    initial_system = build_adaptive_variational_metric_system_v252(
        state, model, settings=settings
    )
    predictor = theta_start + dt * initial_system.velocity

    def residual(theta_end):
        midpoint = 0.5 * (theta_start + np.asarray(theta_end, dtype=float))
        midpoint_state = state_from_adaptive_variational_parameters_v252(
            midpoint,
            ngaussian=state.ngaussian,
            nstate=state.nstate,
            time_au=state.time_au + 0.5 * dt,
        )
        midpoint_system = build_adaptive_variational_metric_system_v252(
            midpoint_state, model, settings=settings
        )
        return (
            np.asarray(theta_end, dtype=float)
            - theta_start
            - dt * midpoint_system.velocity
        )

    predictor_residual_norm = float(np.linalg.norm(residual(predictor)))
    solution = root(
        residual,
        predictor,
        method="hybr",
        options={
            "xtol": float(settings.nonlinear_xtol),
            "maxfev": int(settings.nonlinear_max_function_evaluations),
        },
    )
    theta_end = np.asarray(solution.x, dtype=float)
    final_residual = np.asarray(residual(theta_end), dtype=float)
    final_residual_norm = float(np.linalg.norm(final_residual))
    if not bool(solution.success) or final_residual_norm > settings.nonlinear_residual_tolerance:
        raise RuntimeError(
            "adaptive implicit midpoint TDVP solve failed: "
            f"success={bool(solution.success)}, status={int(solution.status)}, "
            f"residual={final_residual_norm:.6e}, message={solution.message}"
        )
    end = state_from_adaptive_variational_parameters_v252(
        theta_end,
        ngaussian=state.ngaussian,
        nstate=state.nstate,
        time_au=state.time_au + dt,
    ).validate(
        require_normalized=True,
        tolerance=max(
            settings.structural_tolerance, settings.maximum_step_norm_drift
        ),
    )
    _validate_width_domain_v252(end, settings)
    midpoint_parameters = 0.5 * (theta_start + theta_end)
    midpoint_state = state_from_adaptive_variational_parameters_v252(
        midpoint_parameters,
        ngaussian=state.ngaussian,
        nstate=state.nstate,
        time_au=state.time_au + 0.5 * dt,
    )
    midpoint_system = build_adaptive_variational_metric_system_v252(
        midpoint_state, model, settings=settings
    )
    return AdaptiveImplicitMidpointTDVPStepV252(
        start=state,
        end=end,
        dt_au=dt,
        model=model,
        settings=settings,
        midpoint_parameters=midpoint_parameters,
        midpoint_system=midpoint_system,
        nonlinear_success=bool(solution.success),
        nonlinear_status=int(solution.status),
        nonlinear_message=str(solution.message),
        nonlinear_function_evaluations=int(solution.nfev),
        nonlinear_jacobian_evaluations=(
            None if not hasattr(solution, "njev") else int(solution.njev)
        ),
        nonlinear_residual=final_residual,
        nonlinear_residual_norm=final_residual_norm,
        predictor_residual_norm=predictor_residual_norm,
        start_norm=state.generalized_norm,
        end_norm=end.generalized_norm,
        start_energy_hartree=adaptive_variational_energy_v252(state, model),
        end_energy_hartree=adaptive_variational_energy_v252(end, model),
    ).validate()


V252_TDVP_CLAIMS = {
    "adaptive_width_multigaussian_tdvp_validated": True,
    "log_width_positivity_and_quadratic_chirp_validated": True,
    "implicit_midpoint_adaptive_nonlinear_solve_validated": True,
    "svd_metric_rank_and_compatible_nullspace_validated": True,
    "complete_spinor_quadratic_soc_validated": True,
    "single_packet_thawed_harmonic_reduction_validated": True,
    "frozen_coherent_state_reduction_validated": True,
    "gaussian_permutation_covariance_validated": True,
    "constant_electronic_gauge_covariance_validated": True,
    "dynamic_spawning_validated": False,
    "dynamic_pruning_validated": False,
    "coordinate_dependent_electronic_gauge_covariance_validated": False,
    "multidimensional_adaptive_width_tdvp_validated": False,
    "full_correlated_width_matrices_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
    "general_ab_initio_soc_dynamics_accuracy_validated": False,
}


@dataclass(frozen=True)
class AdaptiveWidthMultiGaussianTrajectoryV252:
    initial_state: ThawedGaussianSpinorStateV252
    final_state: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV251
    settings: AdaptiveVariationalSettingsV252
    steps: tuple
    claims: dict = field(default_factory=lambda: dict(V252_TDVP_CLAIMS))

    @property
    def maximum_norm_drift(self):
        initial = self.initial_state.generalized_norm
        return float(max((abs(step.end_norm - initial) for step in self.steps), default=0.0))

    @property
    def maximum_absolute_energy_drift_hartree(self):
        initial = adaptive_variational_energy_v252(self.initial_state, self.model)
        return float(
            max(
                (abs(step.end_energy_hartree - initial) for step in self.steps),
                default=0.0,
            )
        )

    @property
    def minimum_metric_rank(self):
        return int(
            min(
                (step.midpoint_system.solve_receipt.rank for step in self.steps),
                default=self.initial_state.parameter_count,
            )
        )

    @property
    def maximum_metric_nullity(self):
        return int(
            max(
                (step.midpoint_system.solve_receipt.nullity for step in self.steps),
                default=0,
            )
        )

    @property
    def maximum_nonlinear_residual(self):
        return float(
            max(
                (step.nonlinear_residual_norm for step in self.steps), default=0.0
            )
        )

    @property
    def minimum_width(self):
        values = [self.initial_state.widths]
        values.extend(step.end.widths for step in self.steps)
        return float(min(np.min(value) for value in values))

    @property
    def maximum_width(self):
        values = [self.initial_state.widths]
        values.extend(step.end.widths for step in self.steps)
        return float(max(np.max(value) for value in values))

    @property
    def maximum_absolute_chirp(self):
        values = [self.initial_state.chirps]
        values.extend(step.end.chirps for step in self.steps)
        return float(max(np.max(np.abs(value)) for value in values))

    def validate(self):
        model = self.model.validate()
        settings = self.settings.validate()
        tolerance = max(
            settings.structural_tolerance, settings.maximum_step_norm_drift
        )
        self.initial_state.validate(require_normalized=True, tolerance=tolerance)
        self.final_state.validate(require_normalized=True, tolerance=tolerance)
        _validate_width_domain_v252(self.initial_state, settings)
        _validate_width_domain_v252(self.final_state, settings)
        if self.initial_state.nstate != model.nstate:
            raise ValueError("adaptive trajectory state/model dimensions differ.")
        if type(self.claims) is not dict or any(
            type(value) is not bool for value in self.claims.values()
        ):
            raise TypeError("every v0.25.2 TDVP claim must be a native Boolean.")
        if self.claims != V252_TDVP_CLAIMS:
            raise ValueError("v0.25.2 TDVP claims differ from the frozen boundary.")
        previous = self.initial_state
        for step in self.steps:
            step.validate()
            if step.model.fingerprint() != model.fingerprint():
                raise ValueError("adaptive trajectory step model identity changed.")
            if step.settings.as_dict() != settings.as_dict():
                raise ValueError("adaptive trajectory step settings changed.")
            if _scaled_norm_v252(
                pack_adaptive_variational_parameters_v252(step.start),
                pack_adaptive_variational_parameters_v252(previous),
            ) > settings.structural_tolerance or abs(
                step.start.time_au - previous.time_au
            ) > settings.structural_tolerance:
                raise ValueError("v0.25.2 TDVP step chain is discontinuous.")
            previous = step.end
        if _scaled_norm_v252(
            pack_adaptive_variational_parameters_v252(previous),
            pack_adaptive_variational_parameters_v252(self.final_state),
        ) > settings.structural_tolerance or abs(
            previous.time_au - self.final_state.time_au
        ) > settings.structural_tolerance:
            raise ValueError("v0.25.2 TDVP final state is not bound to last step.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": ADAPTIVE_MULTIGAUSSIAN_TDVP_SCHEMA_V252,
            "ansatz": ADAPTIVE_MULTIGAUSSIAN_TDVP_ANSATZ_V252,
            "width_coordinates": WIDTH_COORDINATES_V252,
            "variational_principle": VARIATIONAL_PRINCIPLE_V252,
            "integrator": VARIATIONAL_INTEGRATOR_V252,
            "metric_solver": VARIATIONAL_METRIC_SOLVER_V252,
            "potential_contract": POTENTIAL_CONTRACT_V252,
            "model": self.model.as_dict(),
            "settings": self.settings.as_dict(),
            "initial_state": self.initial_state.as_dict(),
            "final_state": self.final_state.as_dict(),
            "steps": [step.as_dict() for step in self.steps],
            "maximum_norm_drift": self.maximum_norm_drift,
            "maximum_absolute_energy_drift_hartree": (
                self.maximum_absolute_energy_drift_hartree
            ),
            "minimum_metric_rank": self.minimum_metric_rank,
            "maximum_metric_nullity": self.maximum_metric_nullity,
            "maximum_nonlinear_residual": self.maximum_nonlinear_residual,
            "minimum_width": self.minimum_width,
            "maximum_width": self.maximum_width,
            "maximum_absolute_chirp": self.maximum_absolute_chirp,
            "claims": dict(self.claims),
        }

    def fingerprint(self):
        return _sha256_v252(self.as_dict())


def run_adaptive_width_multigaussian_tdvp_v252(
    initial_state,
    model,
    *,
    dt_au,
    steps,
    settings=AdaptiveVariationalSettingsV252(),
):
    settings = settings.validate()
    model = model.validate()
    dt = float(dt_au)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("forward adaptive TDVP dt must be positive.")
    if (
        isinstance(steps, (bool, np.bool_))
        or not isinstance(steps, (int, np.integer))
        or int(steps) < 0
    ):
        raise ValueError("adaptive TDVP steps must be a nonnegative integer.")
    state = initial_state.normalized()
    _validate_width_domain_v252(state, settings)
    initial = state
    receipts = []
    for _ in range(int(steps)):
        receipt = adaptive_implicit_midpoint_tdvp_step_v252(
            state, model, dt, settings=settings
        )
        receipts.append(receipt)
        state = receipt.end
    return AdaptiveWidthMultiGaussianTrajectoryV252(
        initial_state=initial,
        final_state=state,
        model=model,
        settings=settings,
        steps=tuple(receipts),
        claims=dict(V252_TDVP_CLAIMS),
    ).validate()


def reverse_adaptive_width_multigaussian_tdvp_v252(trajectory):
    trajectory = trajectory.validate()
    state = trajectory.final_state
    receipts = []
    for forward_step in reversed(trajectory.steps):
        receipt = adaptive_implicit_midpoint_tdvp_step_v252(
            state,
            trajectory.model,
            -forward_step.dt_au,
            settings=trajectory.settings,
        )
        receipts.append(receipt)
        state = receipt.end
    return AdaptiveWidthMultiGaussianTrajectoryV252(
        initial_state=trajectory.final_state,
        final_state=state,
        model=trajectory.model,
        settings=trajectory.settings,
        steps=tuple(receipts),
        claims=dict(V252_TDVP_CLAIMS),
    ).validate()
