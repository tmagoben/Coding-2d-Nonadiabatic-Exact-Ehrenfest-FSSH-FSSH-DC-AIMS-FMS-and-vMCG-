"""Diagonal-width multidimensional multi-Gaussian TDVP for v0.26.0.

The ansatz is

    Psi_a(R,t) = sum_I C[I,a](t) g_I(R,t),

where every normalized packet has vector centres/momenta and one positive width
and real chirp per coordinate,

    g_I = prod_mu (alpha[I,mu]/pi)**(1/4)
          exp[-(alpha[I,mu]-i beta[I,mu]) (R_mu-q[I,mu])**2/2
              + i p[I,mu] (R_mu-q[I,mu])].

Real parameters are propagated with the McLachlan equations

    G theta_dot = b,
    G_mn = Re <d_m Psi | d_n Psi>,
    b_m  = Im <d_m Psi | H Psi>,

using a full-SVD compatible-null solve and a fully implicit midpoint residual.
All integrals are analytic for the quadratic fixed-frame Hamiltonian contract in
``multidimensional_soc_v260``.  Width matrices remain diagonal in v0.26.0; full
correlated complex width matrices are an explicit future boundary.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np
from scipy.optimize import root

from .multigaussian_tdvp_v251 import MetricSolveReceiptV251, solve_variational_metric_v251
from .multidimensional_soc_v260 import QuadraticSpinHamiltonianNDV260


MULTIDIMENSIONAL_TDVP_SCHEMA_V260 = "gnd-multidimensional-adaptive-gaussian-tdvp-v0.26.0"
MULTIDIMENSIONAL_TDVP_ANSATZ_V260 = (
    "coupled normalized multidimensional Gaussians with coordinate-diagonal logarithmic widths, "
    "coordinate-diagonal chirps, and complete fixed-frame spinors"
)
VARIATIONAL_PRINCIPLE_V260 = "real-parameter McLachlan time-dependent variation"
VARIATIONAL_INTEGRATOR_V260 = "fully implicit midpoint on the SVD-pseudoinverse TDVP vector field"
VARIATIONAL_METRIC_SOLVER_V260 = "full SVD with compatible-null-space audit"
WIDTH_CONVENTION_V260 = "alpha=exp(eta)>0 independently for each packet and coordinate"


def _canonical_tdvp_v260(value):
    if isinstance(value, np.generic):
        return _canonical_tdvp_v260(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_tdvp_v260(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_tdvp_v260(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_tdvp_v260(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("v0.26.0 canonical TDVP data cannot be non-finite.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported TDVP canonical value: {type(value).__name__}")


def _sha256_tdvp_v260(value):
    payload = json.dumps(
        _canonical_tdvp_v260(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_tdvp_v260(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_norm_tdvp_v260(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def _zero_multiindex_v260(ndim):
    return (0,) * int(ndim)


def _monomial_v260(ndim, axes=()):
    powers = [0] * int(ndim)
    for axis in axes:
        powers[int(axis)] += 1
    return tuple(powers)


def _poly_add_v260(*polynomials):
    result = {}
    for polynomial in polynomials:
        for powers, value in polynomial.items():
            result[powers] = result.get(powers, 0.0 + 0.0j) + complex(value)
    return {powers: value for powers, value in result.items() if abs(value) > 0.0}


def _poly_scale_v260(polynomial, scale):
    return {powers: complex(scale) * value for powers, value in polynomial.items()}


def _poly_multiply_v260(left, right):
    result = {}
    for powers_left, value_left in left.items():
        for powers_right, value_right in right.items():
            powers = tuple(a + b for a, b in zip(powers_left, powers_right))
            result[powers] = result.get(powers, 0.0 + 0.0j) + value_left * value_right
    return {powers: value for powers, value in result.items() if abs(value) > 0.0}


def _poly_conjugate_v260(polynomial):
    return {powers: np.conj(value) for powers, value in polynomial.items()}


def _linear_polynomial_v260(ndim, axis, constant, slope):
    return {
        _zero_multiindex_v260(ndim): complex(constant),
        _monomial_v260(ndim, (axis,)): complex(slope),
    }


def _shifted_square_polynomial_v260(ndim, axis, center):
    return {
        _zero_multiindex_v260(ndim): complex(center) ** 2,
        _monomial_v260(ndim, (axis,)): -2.0 * complex(center),
        _monomial_v260(ndim, (axis, axis)): 1.0 + 0.0j,
    }


def _cross_gaussian_data_v260(qi, pi, alphai, betai, qj, pj, alphaj, betaj):
    """Return exact overlap and complex-normal raw-moment parameters."""

    arrays = [np.asarray(value, dtype=float) for value in (qi, pi, alphai, betai, qj, pj, alphaj, betaj)]
    shape = arrays[0].shape
    if len(shape) != 1 or any(value.shape != shape for value in arrays):
        raise ValueError("cross-Gaussian parameters must be equally sized vectors.")
    if not all(np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("cross-Gaussian parameters must be finite.")
    qi, pi, alphai, betai, qj, pj, alphaj, betaj = arrays
    if np.min(alphai) <= 0.0 or np.min(alphaj) <= 0.0:
        raise ValueError("Gaussian widths must be positive.")
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
    overlap = complex(np.exp(np.sum(log_prefactor + constant + 0.5 * linear**2 / combined)))
    mean = linear / combined
    variance = 1.0 / combined
    if not np.isfinite(overlap) or not np.all(np.isfinite(mean)) or not np.all(np.isfinite(variance)):
        raise ValueError("cross-Gaussian moment data are non-finite.")
    return overlap, mean, variance


def _univariate_raw_moments_v260(mean, variance, maximum_order=4):
    values = [1.0 + 0.0j]
    if maximum_order >= 1:
        values.append(mean)
    if maximum_order >= 2:
        values.append(mean**2 + variance)
    if maximum_order >= 3:
        values.append(mean**3 + 3.0 * mean * variance)
    if maximum_order >= 4:
        values.append(mean**4 + 6.0 * mean**2 * variance + 3.0 * variance**2)
    return values


def _moment_table_v260(overlap, mean, variance, maximum_order=4):
    return overlap, tuple(
        tuple(_univariate_raw_moments_v260(mean[axis], variance[axis], maximum_order))
        for axis in range(len(mean))
    )


def _integrate_polynomial_v260(polynomial, moment_table):
    overlap, moments = moment_table
    value = 0.0 + 0.0j
    for powers, coefficient in polynomial.items():
        if len(powers) != len(moments) or max(powers, default=0) >= len(moments[0]):
            raise ValueError("polynomial degree exceeds the available Gaussian moments.")
        raw = overlap
        for axis, power in enumerate(powers):
            raw *= moments[axis][power]
        value += coefficient * raw
    return complex(value)


def _kinetic_polynomial_nd_v260(q, p, widths, chirps, inverse_mass):
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    widths = np.asarray(widths, dtype=float)
    chirps = np.asarray(chirps, dtype=float)
    inverse_mass = np.asarray(inverse_mass, dtype=float)
    ndim = len(q)
    if any(value.shape != (ndim,) for value in (p, widths, chirps)) or inverse_mass.shape != (ndim, ndim):
        raise ValueError("kinetic polynomial dimensions differ.")
    z = widths - 1.0j * chirps
    factors = [
        _linear_polynomial_v260(ndim, axis, z[axis] * q[axis] + 1.0j * p[axis], -z[axis])
        for axis in range(ndim)
    ]
    polynomial = {}
    for a in range(ndim):
        for b in range(ndim):
            contribution = _poly_multiply_v260(factors[a], factors[b])
            if a == b:
                contribution = _poly_add_v260(
                    contribution,
                    {_zero_multiindex_v260(ndim): -z[a]},
                )
            polynomial = _poly_add_v260(
                polynomial,
                _poly_scale_v260(contribution, -0.5 * inverse_mass[a, b]),
            )
    return polynomial


def _operator_polynomials_v260(state, model, packet):
    """Return polynomial matrix entries for H acting on one packet."""

    ndim = state.ndim
    nstate = state.nstate
    result = [[{} for _ in range(nstate)] for _ in range(nstate)]
    zero = _zero_multiindex_v260(ndim)
    for a in range(nstate):
        for b in range(nstate):
            polynomial = {zero: model.H0[a, b]}
            for axis in range(ndim):
                polynomial = _poly_add_v260(
                    polynomial,
                    {_monomial_v260(ndim, (axis,)): model.H1[axis, a, b]},
                )
            for axis_a in range(ndim):
                for axis_b in range(ndim):
                    polynomial = _poly_add_v260(
                        polynomial,
                        {_monomial_v260(ndim, (axis_a, axis_b)): model.H2[axis_a, axis_b, a, b]},
                    )
            result[a][b] = polynomial
    kinetic = _kinetic_polynomial_nd_v260(
        state.q[packet],
        state.p[packet],
        state.widths[packet],
        state.chirps[packet],
        model.inverse_mass_matrix_au,
    )
    for electronic in range(nstate):
        result[electronic][electronic] = _poly_add_v260(result[electronic][electronic], kinetic)
    return result


@dataclass(frozen=True)
class MultidimensionalVariationalSettingsV260:
    metric_relative_cutoff: float = 1.0e-10
    metric_absolute_cutoff: float = 1.0e-12
    maximum_retained_condition_number: float = 1.0e10
    null_rhs_relative_tolerance: float = 3.0e-9
    linear_residual_relative_tolerance: float = 3.0e-9
    nonlinear_residual_tolerance: float = 5.0e-10
    nonlinear_xtol: float = 1.0e-10
    nonlinear_max_function_evaluations: int = 1200
    structural_tolerance: float = 3.0e-10
    maximum_step_norm_drift: float = 8.0e-8
    minimum_width: float = 1.0e-8
    maximum_width: float = 1.0e8
    maximum_absolute_chirp: float = 1.0e8
    maximum_step_log_width_change: float = 0.5
    allow_compatible_rank_deficiency: bool = True
    multidimensional_nuclear_motion: bool = True
    adaptive_diagonal_widths: bool = True
    full_correlated_width_matrices: bool = False
    fixed_electronic_frame: bool = True
    real_molecular_soc_provider: bool = False
    variational_principle: str = VARIATIONAL_PRINCIPLE_V260
    integrator: str = VARIATIONAL_INTEGRATOR_V260
    metric_solver: str = VARIATIONAL_METRIC_SOLVER_V260
    width_coordinates: str = WIDTH_CONVENTION_V260
    nonlinear_solver: str = "scipy.optimize.root-hybr"

    def validate(self):
        positive = (
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
        )
        for name in positive:
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if float(self.metric_relative_cutoff) >= 1.0:
            raise ValueError("metric_relative_cutoff must be smaller than one.")
        if float(self.maximum_retained_condition_number) < 1.0:
            raise ValueError("maximum retained condition number must be at least one.")
        if float(self.maximum_width) <= float(self.minimum_width):
            raise ValueError("maximum_width must exceed minimum_width.")
        if isinstance(self.nonlinear_max_function_evaluations, (bool, np.bool_)) or not isinstance(
            self.nonlinear_max_function_evaluations, (int, np.integer)
        ) or int(self.nonlinear_max_function_evaluations) < 1:
            raise ValueError("nonlinear_max_function_evaluations must be a positive integer.")
        booleans = (
            "allow_compatible_rank_deficiency",
            "multidimensional_nuclear_motion",
            "adaptive_diagonal_widths",
            "full_correlated_width_matrices",
            "fixed_electronic_frame",
            "real_molecular_soc_provider",
        )
        for name in booleans:
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if not self.allow_compatible_rank_deficiency:
            raise ValueError("compatible SVD null-space handling is mandatory.")
        if not self.multidimensional_nuclear_motion or not self.adaptive_diagonal_widths:
            raise ValueError("v0.26.0 freezes multidimensional diagonal adaptive widths as enabled.")
        if self.full_correlated_width_matrices:
            raise ValueError("v0.26.0 does not admit full correlated width matrices.")
        if not self.fixed_electronic_frame:
            raise ValueError("v0.26.0 TDVP requires a fixed electronic frame.")
        if self.real_molecular_soc_provider:
            raise ValueError("v0.26.0 does not admit live molecular-SOC trajectories.")
        if self.variational_principle != VARIATIONAL_PRINCIPLE_V260:
            raise ValueError("the v0.26.0 variational principle is frozen.")
        if self.integrator != VARIATIONAL_INTEGRATOR_V260:
            raise ValueError("the v0.26.0 integrator is frozen.")
        if self.metric_solver != VARIATIONAL_METRIC_SOLVER_V260:
            raise ValueError("the v0.26.0 metric solver is frozen.")
        if self.width_coordinates != WIDTH_CONVENTION_V260:
            raise ValueError("the v0.26.0 width convention is frozen.")
        if self.nonlinear_solver != "scipy.optimize.root-hybr":
            raise ValueError("the v0.26.0 nonlinear solver is frozen.")
        return self

    def as_dict(self):
        self.validate()
        return _canonical_tdvp_v260(asdict(self))


@dataclass(frozen=True)
class DiagonalGaussianSpinorStateV260:
    q: np.ndarray
    p: np.ndarray
    widths: np.ndarray
    chirps: np.ndarray
    coefficients: np.ndarray
    time_au: float = 0.0

    def __post_init__(self):
        for name in ("q", "p", "widths", "chirps"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())
        object.__setattr__(self, "coefficients", np.asarray(self.coefficients, dtype=complex).copy())

    @property
    def ngaussian(self):
        return int(self.q.shape[0]) if self.q.ndim == 2 else 0

    @property
    def ndim(self):
        return int(self.q.shape[1]) if self.q.ndim == 2 else 0

    @property
    def nstate(self):
        return int(self.coefficients.shape[1]) if self.coefficients.ndim == 2 else 0

    @property
    def parameter_count(self):
        return 2 * self.ngaussian * self.nstate + 4 * self.ngaussian * self.ndim

    @property
    def log_widths(self):
        if self.widths.shape != self.q.shape or np.min(self.widths) <= 0.0:
            raise ValueError("cannot form log widths from invalid positive widths.")
        return np.log(self.widths)

    def validate(self, *, require_normalized=False, tolerance=3.0e-10):
        if self.q.ndim != 2 or self.q.shape[0] < 1 or self.q.shape[1] < 1:
            raise ValueError("q must have shape (ngaussian,ndim) with positive dimensions.")
        if any(value.shape != self.q.shape for value in (self.p, self.widths, self.chirps)):
            raise ValueError("all multidimensional packet arrays must share shape.")
        if self.coefficients.ndim != 2 or self.coefficients.shape[0] != self.ngaussian or self.nstate < 1:
            raise ValueError("coefficients must have shape (ngaussian,nstate).")
        if not all(np.all(np.isfinite(value)) for value in (self.q, self.p, self.widths, self.chirps, self.coefficients)):
            raise ValueError("multidimensional Gaussian state contains non-finite values.")
        if np.min(self.widths) <= 0.0:
            raise ValueError("all multidimensional Gaussian widths must be positive.")
        if not np.isfinite(float(self.time_au)):
            raise ValueError("state time must be finite.")
        norm = self.generalized_norm
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError("multidimensional Gaussian state must have positive norm.")
        if require_normalized and abs(norm - 1.0) > float(tolerance):
            raise ValueError("multidimensional Gaussian state is not normalized.")
        return self

    def nuclear_overlap_matrix(self):
        if self.q.ndim != 2 or np.min(self.widths) <= 0.0:
            raise ValueError("cannot build overlap from invalid packet arrays.")
        overlap = np.zeros((self.ngaussian, self.ngaussian), dtype=complex)
        for i in range(self.ngaussian):
            for j in range(self.ngaussian):
                overlap[i, j] = _cross_gaussian_data_v260(
                    self.q[i], self.p[i], self.widths[i], self.chirps[i],
                    self.q[j], self.p[j], self.widths[j], self.chirps[j],
                )[0]
        if _scaled_norm_tdvp_v260(overlap, overlap.conj().T) > 3.0e-10:
            raise ValueError("nuclear overlap matrix is not Hermitian.")
        return overlap

    @property
    def generalized_norm(self):
        if self.coefficients.ndim != 2 or self.q.ndim != 2:
            return float("nan")
        overlap = self.nuclear_overlap_matrix()
        value = np.einsum("ia,ij,ja->", self.coefficients.conj(), overlap, self.coefficients, optimize=True)
        return float(np.real(value))

    def normalized(self):
        self.validate(require_normalized=False)
        return DiagonalGaussianSpinorStateV260(
            self.q, self.p, self.widths, self.chirps,
            self.coefficients / np.sqrt(self.generalized_norm), self.time_au,
        ).validate(require_normalized=True)

    def permuted(self, order):
        self.validate(require_normalized=False)
        order = np.asarray(order, dtype=int)
        if order.shape != (self.ngaussian,) or sorted(order.tolist()) != list(range(self.ngaussian)):
            raise ValueError("packet permutation is invalid.")
        return DiagonalGaussianSpinorStateV260(
            self.q[order], self.p[order], self.widths[order], self.chirps[order],
            self.coefficients[order], self.time_au,
        ).validate(require_normalized=False)

    def gauge_transformed(self, unitary):
        self.validate(require_normalized=False)
        unitary = np.asarray(unitary, dtype=complex)
        if unitary.shape != (self.nstate, self.nstate) or _scaled_norm_tdvp_v260(
            unitary.conj().T @ unitary, np.eye(self.nstate)
        ) > 2.0e-11:
            raise ValueError("electronic gauge transformation must be unitary.")
        return DiagonalGaussianSpinorStateV260(
            self.q,
            self.p,
            self.widths,
            self.chirps,
            self.coefficients @ unitary.conj(),
            self.time_au,
        ).validate(require_normalized=False)

    def coordinate_rotated(self, orthogonal, *, isotropy_tolerance=2.0e-11):
        """Rotate q and p when every diagonal width/chirp block is isotropic.

        A general rotation of an anisotropic diagonal width produces an
        off-diagonal width matrix, which is outside the v0.26.0 manifold and is
        rejected rather than silently discarded.
        """

        self.validate(require_normalized=False)
        orthogonal = np.asarray(orthogonal, dtype=float)
        if orthogonal.shape != (self.ndim, self.ndim) or _scaled_norm_tdvp_v260(
            orthogonal @ orthogonal.T, np.eye(self.ndim)
        ) > float(isotropy_tolerance):
            raise ValueError("coordinate transformation must be orthogonal.")
        signed_permutation = np.allclose(
            np.abs(orthogonal).sum(axis=0),
            1.0,
            atol=isotropy_tolerance,
            rtol=0.0,
        ) and np.allclose(
            np.abs(orthogonal).sum(axis=1),
            1.0,
            atol=isotropy_tolerance,
            rtol=0.0,
        )
        widths = self.widths
        chirps = self.chirps
        if signed_permutation:
            coordinate_permutation = (orthogonal**2).T
            widths = self.widths @ coordinate_permutation
            chirps = self.chirps @ coordinate_permutation
        else:
            if any(
                np.max(np.abs(row - row[0])) > float(isotropy_tolerance)
                for row in self.widths
            ):
                raise ValueError("anisotropic diagonal widths are not closed under general rotations.")
            if any(
                np.max(np.abs(row - row[0])) > float(isotropy_tolerance)
                for row in self.chirps
            ):
                raise ValueError("anisotropic diagonal chirps are not closed under general rotations.")
        return DiagonalGaussianSpinorStateV260(
            self.q @ orthogonal.T,
            self.p @ orthogonal.T,
            widths,
            chirps,
            self.coefficients,
            self.time_au,
        ).validate(require_normalized=False)

    def as_dict(self):
        self.validate(require_normalized=False)
        return {
            "q": self.q.tolist(),
            "p": self.p.tolist(),
            "widths": self.widths.tolist(),
            "chirps": self.chirps.tolist(),
            "coefficients": _complex_pairs_tdvp_v260(self.coefficients),
            "time_au": float(self.time_au),
            "ngaussian": self.ngaussian,
            "ndim": self.ndim,
            "nstate": self.nstate,
            "generalized_norm": self.generalized_norm,
        }


def pack_multidimensional_parameters_v260(state):
    state = state.validate(require_normalized=False)
    return np.concatenate(
        (
            state.coefficients.real.reshape(-1),
            state.coefficients.imag.reshape(-1),
            state.q.reshape(-1),
            state.p.reshape(-1),
            state.log_widths.reshape(-1),
            state.chirps.reshape(-1),
        )
    )


def state_from_multidimensional_parameters_v260(parameters, *, ngaussian, ndim, nstate, time_au):
    parameters = np.asarray(parameters, dtype=float)
    counts = (ngaussian, ndim, nstate)
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) < 1 for value in counts):
        raise ValueError("ngaussian, ndim, and nstate must be positive integers.")
    ngaussian, ndim, nstate = map(int, counts)
    coefficient_count = ngaussian * nstate
    shape_count = ngaussian * ndim
    expected = 2 * coefficient_count + 4 * shape_count
    if parameters.shape != (expected,) or not np.all(np.isfinite(parameters)):
        raise ValueError("multidimensional parameter vector has an invalid shape or values.")
    coefficients = (
        parameters[:coefficient_count] + 1.0j * parameters[coefficient_count : 2 * coefficient_count]
    ).reshape(ngaussian, nstate)
    cursor = 2 * coefficient_count
    q = parameters[cursor : cursor + shape_count].reshape(ngaussian, ndim)
    cursor += shape_count
    p = parameters[cursor : cursor + shape_count].reshape(ngaussian, ndim)
    cursor += shape_count
    eta = parameters[cursor : cursor + shape_count].reshape(ngaussian, ndim)
    cursor += shape_count
    if np.max(np.abs(eta)) > 700.0:
        raise ValueError("a logarithmic width is outside floating-point range.")
    widths = np.exp(eta)
    chirps = parameters[cursor : cursor + shape_count].reshape(ngaussian, ndim)
    return DiagonalGaussianSpinorStateV260(q, p, widths, chirps, coefficients, float(time_au)).validate(
        require_normalized=False
    )


def _validate_width_domain_v260(state, settings):
    state = state.validate(require_normalized=False)
    settings = settings.validate()
    if float(np.min(state.widths)) < float(settings.minimum_width):
        raise ValueError("a Gaussian width fell below the configured minimum.")
    if float(np.max(state.widths)) > float(settings.maximum_width):
        raise ValueError("a Gaussian width exceeded the configured maximum.")
    if float(np.max(np.abs(state.chirps))) > float(settings.maximum_absolute_chirp):
        raise ValueError("a Gaussian chirp exceeded the configured maximum.")
    return state


def build_multidimensional_gaussian_matrices_v260(state, model):
    state = state.validate(require_normalized=False)
    model = model.validate()
    if state.ndim != model.ndim or state.nstate != model.nstate:
        raise ValueError("state dimensions differ from the Hamiltonian model.")
    dimension = state.ngaussian * state.nstate
    overlap = np.zeros((dimension, dimension), dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(state.nstate, dtype=complex)
    for i in range(state.ngaussian):
        si = slice(i * state.nstate, (i + 1) * state.nstate)
        for j in range(state.ngaussian):
            sj = slice(j * state.nstate, (j + 1) * state.nstate)
            data = _cross_gaussian_data_v260(
                state.q[i], state.p[i], state.widths[i], state.chirps[i],
                state.q[j], state.p[j], state.widths[j], state.chirps[j],
            )
            table = _moment_table_v260(*data, maximum_order=2)
            overlap[si, sj] = data[0] * identity
            operator = _operator_polynomials_v260(state, model, j)
            for a in range(state.nstate):
                for b in range(state.nstate):
                    hamiltonian[si.start + a, sj.start + b] = _integrate_polynomial_v260(
                        operator[a][b], table
                    )
    if _scaled_norm_tdvp_v260(overlap, overlap.conj().T) > 3.0e-10:
        raise ValueError("multidimensional overlap matrix is not Hermitian.")
    if _scaled_norm_tdvp_v260(hamiltonian, hamiltonian.conj().T) > 3.0e-10:
        raise ValueError("multidimensional Hamiltonian matrix is not Hermitian.")
    return overlap, hamiltonian


def multidimensional_variational_energy_v260(state, model):
    state = state.validate(require_normalized=False)
    overlap, hamiltonian = build_multidimensional_gaussian_matrices_v260(state, model)
    coefficients = state.coefficients.reshape(-1)
    norm = float(np.real(np.vdot(coefficients, overlap @ coefficients)))
    if norm <= 0.0:
        raise ValueError("cannot evaluate the energy of a zero-norm state.")
    value = np.vdot(coefficients, hamiltonian @ coefficients) / norm
    if abs(float(np.imag(value))) > 3.0e-10:
        raise ValueError("variational energy has a non-negligible imaginary part.")
    return float(np.real(value))


def multidimensional_reduced_density_v260(state, *, normalize=True):
    state = state.validate(require_normalized=False)
    overlap = state.nuclear_overlap_matrix()
    density = state.coefficients.T @ overlap.T @ state.coefficients.conj()
    density = 0.5 * (density + density.conj().T)
    trace = float(np.real(np.trace(density)))
    if normalize:
        if trace <= 0.0:
            raise ValueError("cannot normalize a nonpositive reduced density.")
        density = density / trace
    return density


def _tangent_terms_nd_v260(state):
    """Return (packet, electronic vector, coordinate polynomial) in pack order."""

    state = state.validate(require_normalized=False)
    terms = []
    zero = _zero_multiindex_v260(state.ndim)
    for packet in range(state.ngaussian):
        for electronic in range(state.nstate):
            vector = np.zeros(state.nstate, dtype=complex)
            vector[electronic] = 1.0
            terms.append((packet, vector, {zero: 1.0 + 0.0j}))
    for packet in range(state.ngaussian):
        for electronic in range(state.nstate):
            vector = np.zeros(state.nstate, dtype=complex)
            vector[electronic] = 1.0j
            terms.append((packet, vector, {zero: 1.0 + 0.0j}))
    for packet in range(state.ngaussian):
        z = state.widths[packet] - 1.0j * state.chirps[packet]
        for axis in range(state.ndim):
            terms.append(
                (
                    packet,
                    state.coefficients[packet],
                    _linear_polynomial_v260(
                        state.ndim,
                        axis,
                        -z[axis] * state.q[packet, axis] - 1.0j * state.p[packet, axis],
                        z[axis],
                    ),
                )
            )
    for packet in range(state.ngaussian):
        for axis in range(state.ndim):
            terms.append(
                (
                    packet,
                    state.coefficients[packet],
                    _linear_polynomial_v260(state.ndim, axis, -1.0j * state.q[packet, axis], 1.0j),
                )
            )
    for packet in range(state.ngaussian):
        for axis in range(state.ndim):
            terms.append(
                (
                    packet,
                    state.coefficients[packet],
                    _poly_add_v260(
                        {zero: 0.25 + 0.0j},
                        _poly_scale_v260(
                            _shifted_square_polynomial_v260(state.ndim, axis, state.q[packet, axis]),
                            -0.5 * state.widths[packet, axis],
                        ),
                    ),
                )
            )
    for packet in range(state.ngaussian):
        for axis in range(state.ndim):
            terms.append(
                (
                    packet,
                    state.coefficients[packet],
                    _poly_scale_v260(
                        _shifted_square_polynomial_v260(state.ndim, axis, state.q[packet, axis]), 0.5j
                    ),
                )
            )
    if len(terms) != state.parameter_count:
        raise AssertionError("multidimensional tangent count disagrees with parameter layout.")
    return terms


def active_parameter_indices_v260(state, active_shape_mask):
    state = state.validate(require_normalized=False)
    mask = np.asarray(active_shape_mask, dtype=bool)
    if mask.shape != (state.ngaussian,):
        raise ValueError("active shape mask must contain one value per packet.")
    coefficient_count = 2 * state.ngaussian * state.nstate
    indices = list(range(coefficient_count))
    block = state.ngaussian * state.ndim
    for family in range(4):
        offset = coefficient_count + family * block
        for packet in range(state.ngaussian):
            if mask[packet]:
                indices.extend(range(offset + packet * state.ndim, offset + (packet + 1) * state.ndim))
    return np.asarray(indices, dtype=int)


@dataclass(frozen=True)
class MultidimensionalMetricSystemV260:
    metric: np.ndarray
    rhs: np.ndarray
    velocity: np.ndarray
    active_parameter_indices: np.ndarray
    active_shape_mask: np.ndarray
    solve_receipt: MetricSolveReceiptV251
    generalized_norm: float
    energy_hartree: float
    settings: MultidimensionalVariationalSettingsV260

    def __post_init__(self):
        for name in ("metric", "rhs", "velocity"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())
        object.__setattr__(self, "active_parameter_indices", np.asarray(self.active_parameter_indices, dtype=int).copy())
        object.__setattr__(self, "active_shape_mask", np.asarray(self.active_shape_mask, dtype=bool).copy())

    def validate(self):
        settings = self.settings.validate()
        if self.metric.ndim != 2 or self.metric.shape[0] != self.metric.shape[1]:
            raise ValueError("stored multidimensional metric must be square.")
        active_count = self.metric.shape[0]
        if self.rhs.shape != (active_count,):
            raise ValueError("stored active metric RHS has an incompatible shape.")
        if self.velocity.ndim != 1 or len(self.active_parameter_indices) != active_count:
            raise ValueError("stored full velocity or active index map has an invalid shape.")
        if np.any(self.active_parameter_indices < 0) or np.any(self.active_parameter_indices >= len(self.velocity)):
            raise ValueError("active parameter indices are outside the full vector.")
        if len(np.unique(self.active_parameter_indices)) != active_count:
            raise ValueError("active parameter indices contain duplicates.")
        if not all(np.all(np.isfinite(value)) for value in (self.metric, self.rhs, self.velocity)):
            raise ValueError("stored multidimensional metric system is non-finite.")
        expected_active, expected_receipt = solve_variational_metric_v251(self.metric, self.rhs, settings=settings)
        if _scaled_norm_tdvp_v260(self.velocity[self.active_parameter_indices], expected_active) > settings.structural_tolerance:
            raise ValueError("stored multidimensional velocity disagrees with its SVD solve.")
        inactive = np.ones(len(self.velocity), dtype=bool)
        inactive[self.active_parameter_indices] = False
        if np.max(np.abs(self.velocity[inactive]), initial=0.0) > settings.structural_tolerance:
            raise ValueError("inactive shape coordinates must have zero velocity.")
        if _scaled_norm_tdvp_v260(self.solve_receipt.singular_values, expected_receipt.singular_values) > settings.structural_tolerance:
            raise ValueError("stored multidimensional metric spectrum is inconsistent.")
        if not np.isfinite(float(self.generalized_norm)) or float(self.generalized_norm) <= 0.0:
            raise ValueError("stored generalized norm must be positive.")
        if not np.isfinite(float(self.energy_hartree)):
            raise ValueError("stored variational energy must be finite.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "metric": self.metric.tolist(),
            "rhs": self.rhs.tolist(),
            "velocity": self.velocity.tolist(),
            "active_parameter_indices": self.active_parameter_indices.tolist(),
            "active_shape_mask": self.active_shape_mask.tolist(),
            "solve_receipt": self.solve_receipt.as_dict(),
            "generalized_norm": float(self.generalized_norm),
            "energy_hartree": float(self.energy_hartree),
            "settings": self.settings.as_dict(),
        }


def build_multidimensional_metric_system_v260(
    state,
    model,
    *,
    settings=MultidimensionalVariationalSettingsV260(),
    active_shape_mask=None,
):
    settings = settings.validate()
    state = _validate_width_domain_v260(state, settings)
    model = model.validate()
    if state.ndim != model.ndim or state.nstate != model.nstate:
        raise ValueError("state dimensions differ from the Hamiltonian model.")
    if active_shape_mask is None:
        active_shape_mask = np.ones(state.ngaussian, dtype=bool)
    active_shape_mask = np.asarray(active_shape_mask, dtype=bool)
    active_indices = active_parameter_indices_v260(state, active_shape_mask)
    terms = _tangent_terms_nd_v260(state)
    count = len(terms)
    full_metric = np.zeros((count, count), dtype=float)
    full_rhs = np.zeros(count, dtype=float)
    moment_cache = {}
    operator_cache = {}
    for i in range(state.ngaussian):
        operator_cache[i] = _operator_polynomials_v260(state, model, i)
        for j in range(state.ngaussian):
            data = _cross_gaussian_data_v260(
                state.q[i], state.p[i], state.widths[i], state.chirps[i],
                state.q[j], state.p[j], state.widths[j], state.chirps[j],
            )
            moment_cache[(i, j)] = _moment_table_v260(*data, maximum_order=4)
    for mu, (packet_i, vector_i, polynomial_i) in enumerate(terms):
        conjugate_i = _poly_conjugate_v260(polynomial_i)
        for nu, (packet_j, vector_j, polynomial_j) in enumerate(terms):
            nuclear = _integrate_polynomial_v260(
                _poly_multiply_v260(conjugate_i, polynomial_j), moment_cache[(packet_i, packet_j)]
            )
            full_metric[mu, nu] = float(np.real(np.vdot(vector_i, vector_j) * nuclear))
        projection = 0.0 + 0.0j
        for packet_j in range(state.ngaussian):
            operator = operator_cache[packet_j]
            for output in range(state.nstate):
                electronic_projection = np.conj(vector_i[output])
                if electronic_projection == 0.0:
                    continue
                for input_state in range(state.nstate):
                    amplitude = state.coefficients[packet_j, input_state]
                    if amplitude == 0.0:
                        continue
                    polynomial = _poly_multiply_v260(conjugate_i, operator[output][input_state])
                    projection += (
                        electronic_projection
                        * amplitude
                        * _integrate_polynomial_v260(polynomial, moment_cache[(packet_i, packet_j)])
                    )
        full_rhs[mu] = float(np.imag(projection))
    full_metric = 0.5 * (full_metric + full_metric.T)
    metric = full_metric[np.ix_(active_indices, active_indices)]
    rhs = full_rhs[active_indices]
    active_velocity, receipt = solve_variational_metric_v251(metric, rhs, settings=settings)
    velocity = np.zeros(count, dtype=float)
    velocity[active_indices] = active_velocity
    return MultidimensionalMetricSystemV260(
        metric=metric,
        rhs=rhs,
        velocity=velocity,
        active_parameter_indices=active_indices,
        active_shape_mask=active_shape_mask,
        solve_receipt=receipt,
        generalized_norm=state.generalized_norm,
        energy_hartree=multidimensional_variational_energy_v260(state, model),
        settings=settings,
    ).validate()


@dataclass(frozen=True)
class MultidimensionalImplicitMidpointStepV260:
    start: DiagonalGaussianSpinorStateV260
    end: DiagonalGaussianSpinorStateV260
    model: QuadraticSpinHamiltonianNDV260
    settings: MultidimensionalVariationalSettingsV260
    dt_au: float
    active_shape_mask: np.ndarray
    midpoint_parameters: np.ndarray
    midpoint_system: MultidimensionalMetricSystemV260
    nonlinear_success: bool
    nonlinear_status: int
    nonlinear_message: str
    nonlinear_function_evaluations: int
    nonlinear_residual: np.ndarray
    nonlinear_residual_norm: float
    predictor_residual_norm: float
    start_norm: float
    end_norm: float
    start_energy_hartree: float
    end_energy_hartree: float

    def __post_init__(self):
        object.__setattr__(self, "active_shape_mask", np.asarray(self.active_shape_mask, dtype=bool).copy())
        object.__setattr__(self, "midpoint_parameters", np.asarray(self.midpoint_parameters, dtype=float).copy())
        object.__setattr__(self, "nonlinear_residual", np.asarray(self.nonlinear_residual, dtype=float).copy())

    @property
    def norm_change(self):
        return float(self.end_norm - self.start_norm)

    @property
    def energy_change_hartree(self):
        return float(self.end_energy_hartree - self.start_energy_hartree)

    @property
    def maximum_log_width_change(self):
        return float(np.max(np.abs(np.log(self.end.widths) - np.log(self.start.widths))))

    def validate(self):
        settings = self.settings.validate()
        self.model.validate()
        self.start.validate(require_normalized=True, tolerance=max(settings.maximum_step_norm_drift, settings.structural_tolerance))
        self.end.validate(require_normalized=True, tolerance=max(settings.maximum_step_norm_drift, settings.structural_tolerance))
        if self.start.q.shape != self.end.q.shape or self.start.coefficients.shape != self.end.coefficients.shape:
            raise ValueError("TDVP step changed state dimensions.")
        if self.active_shape_mask.shape != (self.start.ngaussian,):
            raise ValueError("TDVP step active mask has an invalid shape.")
        if not np.isfinite(float(self.dt_au)) or float(self.dt_au) == 0.0:
            raise ValueError("TDVP step dt must be finite and nonzero.")
        if type(self.nonlinear_success) is not bool or not self.nonlinear_success:
            raise ValueError("TDVP nonlinear solve was not successful.")
        if float(self.nonlinear_residual_norm) > float(settings.nonlinear_residual_tolerance):
            raise ValueError("TDVP nonlinear residual exceeds its gate.")
        if abs(self.norm_change) > float(settings.maximum_step_norm_drift):
            raise ValueError("TDVP step norm drift exceeds its gate.")
        if self.maximum_log_width_change > float(settings.maximum_step_log_width_change):
            raise ValueError("TDVP step changed a log width too far.")
        inactive = ~self.active_shape_mask
        if np.any(inactive):
            drift = max(
                float(np.max(np.abs(self.end.q[inactive] - self.start.q[inactive]))),
                float(np.max(np.abs(self.end.p[inactive] - self.start.p[inactive]))),
                float(np.max(np.abs(self.end.widths[inactive] - self.start.widths[inactive]))),
                float(np.max(np.abs(self.end.chirps[inactive] - self.start.chirps[inactive]))),
            )
            if drift > settings.structural_tolerance:
                raise ValueError("inactive newborn packet shapes changed during activation.")
        self.midpoint_system.validate()
        return self

    def as_dict(self):
        self.validate()
        return {
            "start": self.start.as_dict(),
            "end": self.end.as_dict(),
            "dt_au": float(self.dt_au),
            "active_shape_mask": self.active_shape_mask.tolist(),
            "midpoint_parameters": self.midpoint_parameters.tolist(),
            "midpoint_system": self.midpoint_system.as_dict(),
            "nonlinear": {
                "success": self.nonlinear_success,
                "status": int(self.nonlinear_status),
                "message": str(self.nonlinear_message),
                "function_evaluations": int(self.nonlinear_function_evaluations),
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


def multidimensional_implicit_midpoint_step_v260(
    state,
    model,
    dt_au,
    *,
    settings=MultidimensionalVariationalSettingsV260(),
    active_shape_mask=None,
):
    settings = settings.validate()
    model = model.validate()
    state = _validate_width_domain_v260(state, settings).validate(
        require_normalized=True,
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
    )
    if state.ndim != model.ndim or state.nstate != model.nstate:
        raise ValueError("state dimensions differ from the Hamiltonian model.")
    dt = float(dt_au)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError("implicit midpoint dt must be finite and nonzero.")
    if active_shape_mask is None:
        active_shape_mask = np.ones(state.ngaussian, dtype=bool)
    active_shape_mask = np.asarray(active_shape_mask, dtype=bool)
    if active_shape_mask.shape != (state.ngaussian,):
        raise ValueError("active shape mask must contain one value per packet.")
    theta_start = pack_multidimensional_parameters_v260(state)
    initial_system = build_multidimensional_metric_system_v260(
        state, model, settings=settings, active_shape_mask=active_shape_mask
    )
    predictor = theta_start + dt * initial_system.velocity

    def residual(theta_end):
        theta_end = np.asarray(theta_end, dtype=float)
        midpoint = 0.5 * (theta_start + theta_end)
        midpoint_state = state_from_multidimensional_parameters_v260(
            midpoint,
            ngaussian=state.ngaussian,
            ndim=state.ndim,
            nstate=state.nstate,
            time_au=state.time_au + 0.5 * dt,
        )
        midpoint_system = build_multidimensional_metric_system_v260(
            midpoint_state, model, settings=settings, active_shape_mask=active_shape_mask
        )
        return theta_end - theta_start - dt * midpoint_system.velocity

    predictor_residual_norm = float(np.linalg.norm(residual(predictor)))
    solution = root(
        residual,
        predictor,
        method="hybr",
        options={
            "xtol": float(settings.nonlinear_xtol),
            "maxfev": int(settings.nonlinear_max_function_evaluations),
            # Bound MINPACK's initial trust-region radius.  Multidimensional
            # log-width coordinates are mathematically unbounded, so the
            # default factor=100 can probe meaningless exp(eta) values before
            # the already accurate explicit predictor is accepted.
            "factor": 0.1,
            "diag": np.maximum(np.abs(theta_start), 1.0),
        },
    )
    theta_end = np.asarray(solution.x, dtype=float)
    # MINPACK solves in the full real parameter vector.  Inactive shape
    # equations are mathematically ``theta_end - theta_start = 0``, but the
    # nonlinear iteration can leave round-off-sized values in those slots.
    # Restore the frozen coordinates explicitly so coefficient-only newborn
    # propagation is bitwise contractual, not merely tolerance-based.
    active_indices = active_parameter_indices_v260(state, active_shape_mask)
    frozen = np.ones(theta_end.size, dtype=bool)
    frozen[active_indices] = False
    theta_end[frozen] = theta_start[frozen]
    nonlinear_residual = np.asarray(residual(theta_end), dtype=float)
    nonlinear_residual_norm = float(np.linalg.norm(nonlinear_residual))
    if not bool(solution.success) or nonlinear_residual_norm > settings.nonlinear_residual_tolerance:
        raise RuntimeError(
            "multidimensional implicit midpoint TDVP solve failed: "
            f"success={bool(solution.success)}, status={int(solution.status)}, "
            f"residual={nonlinear_residual_norm:.6e}, message={solution.message}"
        )
    end = state_from_multidimensional_parameters_v260(
        theta_end,
        ngaussian=state.ngaussian,
        ndim=state.ndim,
        nstate=state.nstate,
        time_au=state.time_au + dt,
    ).validate(
        require_normalized=True,
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
    )
    _validate_width_domain_v260(end, settings)
    midpoint_parameters = 0.5 * (theta_start + theta_end)
    midpoint_state = state_from_multidimensional_parameters_v260(
        midpoint_parameters,
        ngaussian=state.ngaussian,
        ndim=state.ndim,
        nstate=state.nstate,
        time_au=state.time_au + 0.5 * dt,
    )
    midpoint_system = build_multidimensional_metric_system_v260(
        midpoint_state, model, settings=settings, active_shape_mask=active_shape_mask
    )
    return MultidimensionalImplicitMidpointStepV260(
        start=state,
        end=end,
        model=model,
        settings=settings,
        dt_au=dt,
        active_shape_mask=active_shape_mask,
        midpoint_parameters=midpoint_parameters,
        midpoint_system=midpoint_system,
        nonlinear_success=bool(solution.success),
        nonlinear_status=int(solution.status),
        nonlinear_message=str(solution.message),
        nonlinear_function_evaluations=int(solution.nfev),
        nonlinear_residual=nonlinear_residual,
        nonlinear_residual_norm=nonlinear_residual_norm,
        predictor_residual_norm=predictor_residual_norm,
        start_norm=state.generalized_norm,
        end_norm=end.generalized_norm,
        start_energy_hartree=multidimensional_variational_energy_v260(state, model),
        end_energy_hartree=multidimensional_variational_energy_v260(end, model),
    ).validate()


@dataclass(frozen=True)
class MultidimensionalTDVPTrajectoryV260:
    initial_state: DiagonalGaussianSpinorStateV260
    final_state: DiagonalGaussianSpinorStateV260
    model: QuadraticSpinHamiltonianNDV260
    settings: MultidimensionalVariationalSettingsV260
    steps: tuple

    @property
    def maximum_norm_drift(self):
        values = [abs(step.norm_change) for step in self.steps]
        values.append(abs(self.final_state.generalized_norm - self.initial_state.generalized_norm))
        return float(max(values, default=0.0))

    @property
    def maximum_energy_drift_hartree(self):
        initial = multidimensional_variational_energy_v260(self.initial_state, self.model)
        values = [abs(multidimensional_variational_energy_v260(step.end, self.model) - initial) for step in self.steps]
        return float(max(values, default=0.0))

    def validate(self):
        self.model.validate()
        settings = self.settings.validate()
        self.initial_state.validate(require_normalized=True, tolerance=settings.maximum_step_norm_drift)
        self.final_state.validate(require_normalized=True, tolerance=settings.maximum_step_norm_drift)
        previous = self.initial_state
        for step in self.steps:
            step.validate()
            if _scaled_norm_tdvp_v260(pack_multidimensional_parameters_v260(previous), pack_multidimensional_parameters_v260(step.start)) > settings.structural_tolerance:
                raise ValueError("TDVP trajectory steps are not contiguous.")
            previous = step.end
        if _scaled_norm_tdvp_v260(pack_multidimensional_parameters_v260(previous), pack_multidimensional_parameters_v260(self.final_state)) > settings.structural_tolerance:
            raise ValueError("TDVP final state differs from the final step.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIDIMENSIONAL_TDVP_SCHEMA_V260,
            "ansatz": MULTIDIMENSIONAL_TDVP_ANSATZ_V260,
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "initial_state": self.initial_state.as_dict(),
            "final_state": self.final_state.as_dict(),
            "steps": [step.as_dict() for step in self.steps],
            "maximum_norm_drift": self.maximum_norm_drift,
            "maximum_energy_drift_hartree": self.maximum_energy_drift_hartree,
        }

    def fingerprint(self):
        return _sha256_tdvp_v260(self.as_dict())


def run_multidimensional_tdvp_v260(
    initial_state,
    model,
    dt_au,
    steps,
    *,
    settings=MultidimensionalVariationalSettingsV260(),
    active_shape_mask=None,
):
    settings = settings.validate()
    initial_state = initial_state.validate(require_normalized=True, tolerance=settings.maximum_step_norm_drift)
    if isinstance(steps, (bool, np.bool_)) or not isinstance(steps, (int, np.integer)) or int(steps) < 0:
        raise ValueError("steps must be a nonnegative integer.")
    current = initial_state
    receipts = []
    for _ in range(int(steps)):
        receipt = multidimensional_implicit_midpoint_step_v260(
            current,
            model,
            dt_au,
            settings=settings,
            active_shape_mask=active_shape_mask,
        )
        receipts.append(receipt)
        current = receipt.end
    return MultidimensionalTDVPTrajectoryV260(
        initial_state=initial_state,
        final_state=current,
        model=model,
        settings=settings,
        steps=tuple(receipts),
    ).validate()


def evaluate_multidimensional_state_v260(state, points):
    """Evaluate Psi at points[...,ndim], returning points_shape+(nstate,)."""

    state = state.validate(require_normalized=False)
    points = np.asarray(points, dtype=float)
    if points.shape == () or points.shape[-1] != state.ndim or not np.all(np.isfinite(points)):
        raise ValueError("evaluation points must be finite with final axis ndim.")
    wavefunction = np.zeros(points.shape[:-1] + (state.nstate,), dtype=complex)
    for packet in range(state.ngaussian):
        displacement = points - state.q[packet]
        normalization = float(np.prod(state.widths[packet] / np.pi) ** 0.25)
        exponent = np.sum(
            -0.5 * state.widths[packet] * displacement**2
            + 0.5j * state.chirps[packet] * displacement**2
            + 1.0j * state.p[packet] * displacement,
            axis=-1,
        )
        gaussian = normalization * np.exp(exponent)
        wavefunction += gaussian[..., None] * state.coefficients[packet]
    return wavefunction


def multidimensional_state_on_grid_v260(state, grid):
    state = state.validate(require_normalized=False)
    points = grid.validate().mesh()
    values = evaluate_multidimensional_state_v260(state, points)
    return np.moveaxis(values, -1, 0)


def residual_coupling_at_geometry_v260(
    state,
    model,
    velocity,
    *,
    q,
    p,
    widths,
    chirps,
):
    """Return ``<g_candidate,a|dPsi/dt+i H Psi>`` for every electronic state.

    This coefficient-tangent residual is used by the controlled basis layer.  It
    is evaluated analytically and does not sample a grid.
    """

    state = state.validate(require_normalized=False)
    model = model.validate()
    velocity = np.asarray(velocity, dtype=float)
    if velocity.shape != (state.parameter_count,) or not np.all(np.isfinite(velocity)):
        raise ValueError("TDVP velocity has an invalid shape or values.")
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    widths = np.asarray(widths, dtype=float)
    chirps = np.asarray(chirps, dtype=float)
    if any(value.shape != (state.ndim,) for value in (q, p, widths, chirps)):
        raise ValueError("candidate geometry vectors must match state dimensionality.")
    if not all(np.all(np.isfinite(value)) for value in (q, p, widths, chirps)) or np.min(widths) <= 0.0:
        raise ValueError("candidate geometry must be finite with positive widths.")
    tables = {}
    for packet in range(state.ngaussian):
        data = _cross_gaussian_data_v260(
            q,
            p,
            widths,
            chirps,
            state.q[packet],
            state.p[packet],
            state.widths[packet],
            state.chirps[packet],
        )
        tables[packet] = _moment_table_v260(*data, maximum_order=4)
    derivative = np.zeros(state.nstate, dtype=complex)
    for speed, (packet, electronic_vector, polynomial) in zip(velocity, _tangent_terms_nd_v260(state)):
        if speed != 0.0:
            derivative += float(speed) * electronic_vector * _integrate_polynomial_v260(
                polynomial, tables[packet]
            )
    hamiltonian_action = np.zeros(state.nstate, dtype=complex)
    for packet in range(state.ngaussian):
        operator = _operator_polynomials_v260(state, model, packet)
        for output in range(state.nstate):
            for input_state in range(state.nstate):
                amplitude = state.coefficients[packet, input_state]
                if amplitude != 0.0:
                    hamiltonian_action[output] += amplitude * _integrate_polynomial_v260(
                        operator[output][input_state], tables[packet]
                    )
    return derivative + 1.0j * hamiltonian_action


V260_MULTIDIMENSIONAL_TDVP_CLAIMS = {
    "multidimensional_adaptive_width_tdvp_validated": True,
    "coordinate_diagonal_log_widths_and_chirps_validated": True,
    "full_positive_definite_mass_matrix_validated": True,
    "implicit_midpoint_nonlinear_solve_validated": True,
    "svd_metric_null_space_handling_validated": True,
    "complete_spinor_complex_soc_validated": True,
    "packet_permutation_covariance_validated": True,
    "constant_electronic_gauge_covariance_validated": True,
    "one_dimensional_v0252_reduction_validated": True,
    "full_correlated_width_matrices_validated": False,
    "coordinate_dependent_electronic_gauge_covariance_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
    "general_ab_initio_soc_accuracy_validated": False,
}
