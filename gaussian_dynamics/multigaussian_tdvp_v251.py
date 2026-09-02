"""Frozen-width multi-Gaussian TDVP with an implicit midpoint solver.

v0.25.1 advances the v0.25.0 single-packet restriction to a genuinely coupled
multi-Gaussian variational metric.  The released validation contract is deliberately
narrow and exact: one nuclear coordinate, fixed Gaussian widths, a complete
electronic spinor in one fixed frame, and a Hermitian matrix potential that is at
most quadratic in the coordinate.

For real parameters ``theta`` the McLachlan equations are

    G_{mu nu} theta_dot_nu = b_mu,
    G_{mu nu} = Re <d_mu Psi | d_nu Psi>,
    b_mu = Im <d_mu Psi | H Psi>.

All Gaussian overlap, tangent-metric, kinetic, and quadratic-potential moments are
analytic.  Rank and compatible null directions are handled by a full SVD.  Time
propagation solves the nonlinear implicit-midpoint residual; it does not apply
independent velocity-Verlet updates to the Gaussian centers.

Adaptive widths, spawning, pruning, coordinate-dependent electronic frames,
multidimensional nuclear motion, and real molecular-SOC trajectories remain closed.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np
from scipy.optimize import root

from .gaussian_general import (
    gaussian_cross_centroid,
    gaussian_cross_covariance,
    gaussian_overlap_general,
)


MULTIGAUSSIAN_TDVP_SCHEMA_V251 = (
    "gnd-frozen-width-multigaussian-tdvp-trajectory-v0.25.1"
)
MULTIGAUSSIAN_TDVP_ANSATZ_V251 = (
    "one-dimensional fixed-width coupled multi-Gaussian packets with complete spinors"
)
VARIATIONAL_PRINCIPLE_V251 = "real-parameter McLachlan time-dependent variation"
VARIATIONAL_INTEGRATOR_V251 = (
    "fully implicit midpoint applied to the SVD-pseudoinverse TDVP vector field"
)
VARIATIONAL_METRIC_SOLVER_V251 = (
    "full SVD with relative/absolute rank cutoff and compatible-null-space audit"
)
POTENTIAL_CONTRACT_V251 = (
    "fixed-frame one-dimensional Hermitian quadratic complete-spinor Hamiltonian"
)


def _canonical_v251(value):
    if isinstance(value, np.generic):
        return _canonical_v251(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_v251(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_v251(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v251(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("v0.25.1 canonical data cannot contain non-finite values.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.25.1 canonical value: {type(value).__name__}")


def _sha256_v251(value):
    payload = json.dumps(
        _canonical_v251(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_v251(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_norm_v251(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def _hermiticity_residual_v251(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    return _scaled_norm_v251(matrix, matrix.conj().T)


def _unitary_v251(matrix, *, tolerance):
    matrix = np.asarray(matrix, dtype=complex)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("electronic gauge matrix must be nonempty and square.")
    identity = np.eye(matrix.shape[0], dtype=complex)
    if _scaled_norm_v251(matrix.conj().T @ matrix, identity) > tolerance:
        raise ValueError("electronic gauge matrix must be unitary.")
    return matrix


@dataclass(frozen=True)
class VariationalMetricSettingsV251:
    """Frozen algorithms and numerical gates for the v0.25.1 metric solve."""

    metric_relative_cutoff: float = 1.0e-10
    metric_absolute_cutoff: float = 1.0e-12
    maximum_retained_condition_number: float = 1.0e10
    null_rhs_relative_tolerance: float = 2.0e-9
    linear_residual_relative_tolerance: float = 2.0e-9
    nonlinear_residual_tolerance: float = 2.0e-10
    nonlinear_xtol: float = 1.0e-11
    nonlinear_max_function_evaluations: int = 600
    structural_tolerance: float = 2.0e-10
    maximum_step_norm_drift: float = 2.0e-8
    variational_principle: str = VARIATIONAL_PRINCIPLE_V251
    integrator: str = VARIATIONAL_INTEGRATOR_V251
    metric_solver: str = VARIATIONAL_METRIC_SOLVER_V251
    nonlinear_solver: str = "scipy.optimize.root-hybr"
    allow_compatible_rank_deficiency: bool = True
    adaptive_gaussian_widths: bool = False
    spawning: bool = False
    pruning: bool = False
    coordinate_dependent_electronic_frame: bool = False
    multidimensional_nuclear_motion: bool = False
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
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if float(self.metric_relative_cutoff) >= 1.0:
            raise ValueError("metric_relative_cutoff must be smaller than one.")
        if float(self.maximum_retained_condition_number) < 1.0:
            raise ValueError("maximum_retained_condition_number must be at least one.")
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
            "real_molecular_soc_provider",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if not self.allow_compatible_rank_deficiency:
            raise ValueError(
                "v0.25.1 freezes compatible SVD null-space handling as enabled."
            )
        closed = {
            "adaptive_gaussian_widths": self.adaptive_gaussian_widths,
            "spawning": self.spawning,
            "pruning": self.pruning,
            "coordinate_dependent_electronic_frame": (
                self.coordinate_dependent_electronic_frame
            ),
            "multidimensional_nuclear_motion": self.multidimensional_nuclear_motion,
            "real_molecular_soc_provider": self.real_molecular_soc_provider,
        }
        requested = [name for name, enabled in closed.items() if enabled]
        if requested:
            raise ValueError(
                "v0.25.1 does not admit: " + ", ".join(sorted(requested)) + "."
            )
        if self.variational_principle != VARIATIONAL_PRINCIPLE_V251:
            raise ValueError("the v0.25.1 variational principle is frozen.")
        if self.integrator != VARIATIONAL_INTEGRATOR_V251:
            raise ValueError("the v0.25.1 implicit integrator is frozen.")
        if self.metric_solver != VARIATIONAL_METRIC_SOLVER_V251:
            raise ValueError("the v0.25.1 SVD metric solver is frozen.")
        if self.nonlinear_solver != "scipy.optimize.root-hybr":
            raise ValueError("the v0.25.1 nonlinear solver is frozen.")
        return self

    def as_dict(self):
        self.validate()
        return _canonical_v251(asdict(self))


@dataclass(frozen=True)
class QuadraticSpinHamiltonianV251:
    """One-dimensional fixed-frame matrix Hamiltonian H(x)=H0+x H1+x^2 H2."""

    mass_au: float
    H0: np.ndarray
    H1: np.ndarray
    H2: np.ndarray
    label: str = "quadratic complete-spinor Hamiltonian"
    physical_soc: bool = False
    complete_spin_manifold: bool = True
    source: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ("H0", "H1", "H2"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=complex).copy()
            )
        object.__setattr__(self, "source", dict(self.source))

    @property
    def nstate(self):
        return int(self.H0.shape[0]) if self.H0.ndim == 2 else 0

    def validate(self, tolerance=2.0e-11):
        mass = float(self.mass_au)
        if not np.isfinite(mass) or mass <= 0.0:
            raise ValueError("quadratic spin Hamiltonian mass must be positive.")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("quadratic spin Hamiltonian label must be nonempty.")
        for name in ("physical_soc", "complete_spin_manifold"):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if not self.complete_spin_manifold:
            raise ValueError("v0.25.1 requires a complete electronic spin manifold.")
        if self.H0.ndim != 2 or self.H0.shape[0] < 1 or self.H0.shape[0] != self.H0.shape[1]:
            raise ValueError("H0 must be a nonempty square matrix.")
        if self.H1.shape != self.H0.shape or self.H2.shape != self.H0.shape:
            raise ValueError("quadratic Hamiltonian coefficient shapes differ.")
        if any(
            not np.all(np.isfinite(matrix))
            for matrix in (self.H0, self.H1, self.H2)
        ):
            raise ValueError("quadratic Hamiltonian contains non-finite values.")
        for name, matrix in (("H0", self.H0), ("H1", self.H1), ("H2", self.H2)):
            if _hermiticity_residual_v251(matrix) > float(tolerance):
                raise ValueError(f"{name} must be Hermitian.")
        _sha256_v251(self.source)
        return self

    def hamiltonian(self, coordinate):
        self.validate()
        x = float(coordinate)
        if not np.isfinite(x):
            raise ValueError("quadratic Hamiltonian coordinate must be finite.")
        return self.H0 + x * self.H1 + x * x * self.H2

    def derivative(self, coordinate):
        self.validate()
        x = float(coordinate)
        if not np.isfinite(x):
            raise ValueError("quadratic Hamiltonian coordinate must be finite.")
        return self.H1 + 2.0 * x * self.H2

    def gauge_transformed(self, unitary):
        unitary = _unitary_v251(unitary, tolerance=2.0e-11)
        if unitary.shape != (self.nstate, self.nstate):
            raise ValueError("electronic gauge dimension disagrees with the model.")
        return QuadraticSpinHamiltonianV251(
            mass_au=self.mass_au,
            H0=unitary.conj().T @ self.H0 @ unitary,
            H1=unitary.conj().T @ self.H1 @ unitary,
            H2=unitary.conj().T @ self.H2 @ unitary,
            label=self.label + " [constant gauge]",
            physical_soc=self.physical_soc,
            complete_spin_manifold=self.complete_spin_manifold,
            source={
                **self.source,
                "constant_gauge_transform": _complex_pairs_v251(unitary),
            },
        ).validate()

    def as_dict(self):
        self.validate()
        return {
            "contract": POTENTIAL_CONTRACT_V251,
            "mass_au": float(self.mass_au),
            "H0": _complex_pairs_v251(self.H0),
            "H1": _complex_pairs_v251(self.H1),
            "H2": _complex_pairs_v251(self.H2),
            "label": self.label,
            "physical_soc": self.physical_soc,
            "complete_spin_manifold": self.complete_spin_manifold,
            "source": _canonical_v251(self.source),
        }

    def fingerprint(self):
        return _sha256_v251(self.as_dict())


def quadratic_spin_hamiltonian_from_provider_v251(
    provider,
    *,
    label=None,
    fitting_coordinates=(-1.0, 0.0, 1.0),
    audit_coordinates=(-1.4, -0.37, 0.29, 1.23),
    tolerance=2.0e-10,
):
    """Freeze a verified quadratic fixed-frame model from an analytic provider.

    The provider must expose complete arbitrary-geometry snapshots, a constant mass,
    zero connection, and one common electronic frame.  Static/differential PySCF
    receipts therefore fail this contract.
    """

    if not callable(getattr(provider, "evaluate_snapshot", None)):
        raise TypeError("quadratic TDVP intake requires evaluate_snapshot(q).")
    provenance = getattr(provider, "provenance", None)
    if provenance is None or not callable(getattr(provenance, "validate", None)):
        raise TypeError("quadratic TDVP intake requires explicit operator provenance.")
    provenance = provenance.validate()
    if not provenance.model_space.fixed_frame:
        raise ValueError("quadratic TDVP intake requires fixed-frame provenance.")
    if not provenance.model_space.complete_multiplets:
        raise ValueError("quadratic TDVP intake requires complete spin multiplets.")
    fitting_coordinates = tuple(float(x) for x in fitting_coordinates)
    if fitting_coordinates != (-1.0, 0.0, 1.0):
        raise ValueError("v0.25.1 freezes fitting coordinates to (-1,0,1).")
    snapshots = [
        provider.evaluate_snapshot(np.asarray([x], dtype=float)).validate()
        for x in fitting_coordinates
    ]
    points = [snapshot.point for snapshot in snapshots]
    nstate = points[0].nstate
    reference_vectors = np.asarray(snapshots[0].state_vectors, dtype=complex)
    masses = []
    for snapshot, point in zip(snapshots, points):
        if point.nq != 1 or point.nstate != nstate:
            raise ValueError("v0.25.1 quadratic intake requires one coordinate and fixed state dimension.")
        if point.mass_matrix_q_au.shape != (1, 1):
            raise ValueError("quadratic TDVP intake requires a scalar mass.")
        masses.append(float(point.mass_matrix_q_au[0, 0]))
        if np.max(np.abs(point.connection_q)) > float(tolerance):
            raise ValueError("v0.25.1 quadratic TDVP intake requires a fixed electronic frame.")
        if not np.allclose(
            snapshot.state_vectors,
            reference_vectors,
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("provider electronic frame changes across fitting coordinates.")
    if not np.allclose(masses, masses[0], atol=tolerance, rtol=tolerance):
        raise ValueError("quadratic TDVP intake requires constant mass.")
    H_minus, H_zero, H_plus = (np.asarray(point.H, dtype=complex) for point in points)
    H0 = H_zero
    H1 = 0.5 * (H_plus - H_minus)
    H2 = 0.5 * (H_plus + H_minus) - H_zero
    metadata = dict(getattr(snapshots[0], "metadata", {}))
    model = QuadraticSpinHamiltonianV251(
        mass_au=masses[0],
        H0=H0,
        H1=H1,
        H2=H2,
        label=(
            str(label)
            if label is not None
            else str(metadata.get("provider", type(provider).__name__))
        ),
        physical_soc=bool(metadata.get("physical_soc", False)),
        complete_spin_manifold=bool(provenance.model_space.complete_multiplets),
        source={
            "provider_type": type(provider).__name__,
            "fitting_coordinates": list(fitting_coordinates),
            "provider_provenance_fingerprint": metadata.get(
                "provenance_fingerprint"
            ),
            "model_space": provenance.model_space.as_dict(),
            "fixed_frame_verified": True,
            "quadratic_H_and_linear_K_verified": True,
        },
    ).validate(tolerance=tolerance)
    for x in tuple(float(item) for item in audit_coordinates):
        snapshot = provider.evaluate_snapshot(np.asarray([x], dtype=float)).validate()
        point = snapshot.point
        if not np.allclose(
            point.mass_matrix_q_au,
            np.asarray([[model.mass_au]]),
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("provider mass changes on the quadratic audit stencil.")
        if np.max(np.abs(point.connection_q)) > float(tolerance):
            raise ValueError("provider connection is nonzero on the quadratic audit stencil.")
        if not np.allclose(
            snapshot.state_vectors,
            reference_vectors,
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("provider frame changes on the quadratic audit stencil.")
        if not np.allclose(
            point.H,
            model.hamiltonian(x),
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("provider Hamiltonian is not quadratic on the audit stencil.")
        if not np.allclose(
            point.hamiltonian_derivative_operator_q[0],
            model.derivative(x),
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("provider physical derivative is not the quadratic-model derivative.")
    return model


def _cross_moments_v251(qi, pi, ai, qj, pj, aj, maximum_order=3):
    if int(maximum_order) != maximum_order or not 0 <= int(maximum_order) <= 3:
        raise ValueError("v0.25.1 cross moments are implemented through order three.")
    Ai = np.asarray([[float(ai)]])
    Aj = np.asarray([[float(aj)]])
    overlap = gaussian_overlap_general(
        np.asarray([qi]),
        np.asarray([pi]),
        Ai,
        np.asarray([qj]),
        np.asarray([pj]),
        Aj,
    )
    mean = gaussian_cross_centroid(
        np.asarray([qi]),
        np.asarray([pi]),
        Ai,
        np.asarray([qj]),
        np.asarray([pj]),
        Aj,
    )[0]
    variance = gaussian_cross_covariance(Ai, Aj)[0, 0]
    values = [overlap]
    if maximum_order >= 1:
        values.append(overlap * mean)
    if maximum_order >= 2:
        values.append(overlap * (mean * mean + variance))
    if maximum_order >= 3:
        values.append(overlap * (mean**3 + 3.0 * mean * variance))
    return np.asarray(values, dtype=complex)


def _integrate_polynomial_v251(coefficients, moments):
    coefficients = np.asarray(coefficients, dtype=complex)
    moments = np.asarray(moments, dtype=complex)
    if coefficients.ndim != 1 or len(coefficients) > len(moments):
        raise ValueError("polynomial degree exceeds the available Gaussian moments.")
    return complex(np.dot(coefficients, moments[: len(coefficients)]))


def _kinetic_polynomial_v251(q, p, width, mass):
    q = float(q)
    p = float(p)
    width = float(width)
    mass = float(mass)
    return np.asarray(
        [
            (-width**2 * q**2 - 2.0j * width * p * q + p**2 + width)
            / (2.0 * mass),
            (width**2 * q + 1.0j * width * p) / mass,
            -(width**2) / (2.0 * mass),
        ],
        dtype=complex,
    )


@dataclass(frozen=True)
class FrozenGaussianSpinorStateV251:
    q: np.ndarray
    p: np.ndarray
    widths: np.ndarray
    coefficients: np.ndarray
    time_au: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "q", np.asarray(self.q, dtype=float).copy())
        object.__setattr__(self, "p", np.asarray(self.p, dtype=float).copy())
        object.__setattr__(
            self, "widths", np.asarray(self.widths, dtype=float).copy()
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
        return 2 * self.ngaussian * self.nstate + 2 * self.ngaussian

    def nuclear_overlap_matrix(self):
        if (
            self.q.ndim != 1
            or len(self.q) < 1
            or self.p.shape != self.q.shape
            or self.widths.shape != self.q.shape
            or np.any(self.widths <= 0.0)
        ):
            raise ValueError("cannot build overlap from an invalid Gaussian geometry.")
        overlap = np.zeros((self.ngaussian, self.ngaussian), dtype=complex)
        for i in range(self.ngaussian):
            for j in range(self.ngaussian):
                overlap[i, j] = _cross_moments_v251(
                    self.q[i],
                    self.p[i],
                    self.widths[i],
                    self.q[j],
                    self.p[j],
                    self.widths[j],
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
            raise ValueError("multi-Gaussian coordinates must be a nonempty vector.")
        if self.p.shape != self.q.shape or self.widths.shape != self.q.shape:
            raise ValueError("multi-Gaussian q, p, and width shapes differ.")
        if self.coefficients.ndim != 2 or self.coefficients.shape[0] != len(self.q):
            raise ValueError("coefficient array must have shape (n_gaussian,n_state).")
        if self.coefficients.shape[1] < 1:
            raise ValueError("multi-Gaussian state requires at least one electronic state.")
        if not all(
            np.all(np.isfinite(value))
            for value in (self.q, self.p, self.widths, self.coefficients)
        ):
            raise ValueError("multi-Gaussian state contains non-finite values.")
        if np.any(self.widths <= 0.0):
            raise ValueError("all frozen Gaussian widths must be positive.")
        if not np.isfinite(float(self.time_au)):
            raise ValueError("multi-Gaussian time must be finite.")
        norm = self.generalized_norm
        if not np.isfinite(norm) or norm <= tolerance:
            raise ValueError("multi-Gaussian wavefunction norm must be positive.")
        if require_normalized and abs(norm - 1.0) > tolerance:
            raise ValueError("multi-Gaussian wavefunction must have unit generalized norm.")
        return self

    def normalized(self):
        self.validate(require_normalized=False)
        return FrozenGaussianSpinorStateV251(
            q=self.q,
            p=self.p,
            widths=self.widths,
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
            raise ValueError("Gaussian permutation must contain every packet exactly once.")
        return FrozenGaussianSpinorStateV251(
            q=self.q[order],
            p=self.p[order],
            widths=self.widths[order],
            coefficients=self.coefficients[order],
            time_au=self.time_au,
        ).validate(tolerance=2.0e-9)

    def gauge_transformed(self, unitary):
        unitary = _unitary_v251(unitary, tolerance=2.0e-11)
        if unitary.shape != (self.nstate, self.nstate):
            raise ValueError("electronic gauge dimension disagrees with the state.")
        return FrozenGaussianSpinorStateV251(
            q=self.q,
            p=self.p,
            widths=self.widths,
            coefficients=self.coefficients @ unitary.conj(),
            time_au=self.time_au,
        ).validate(tolerance=2.0e-9)

    def as_dict(self):
        self.validate(require_normalized=False)
        return {
            "q": self.q.tolist(),
            "p": self.p.tolist(),
            "widths": self.widths.tolist(),
            "coefficients": _complex_pairs_v251(self.coefficients),
            "time_au": float(self.time_au),
            "generalized_norm": self.generalized_norm,
            "n_gaussian": self.ngaussian,
            "n_state": self.nstate,
            "parameter_count": self.parameter_count,
        }


def pack_variational_parameters_v251(state):
    state = state.validate(require_normalized=False)
    return np.concatenate(
        (
            state.coefficients.real.reshape(-1),
            state.coefficients.imag.reshape(-1),
            state.q,
            state.p,
        )
    )


def state_from_variational_parameters_v251(
    parameters,
    *,
    widths,
    nstate,
    time_au,
):
    parameters = np.asarray(parameters, dtype=float)
    widths = np.asarray(widths, dtype=float)
    nstate = int(nstate)
    ngaussian = len(widths)
    expected = 2 * ngaussian * nstate + 2 * ngaussian
    if parameters.shape != (expected,):
        raise ValueError("variational parameter vector has incompatible length.")
    block = ngaussian * nstate
    coefficients = (
        parameters[:block] + 1.0j * parameters[block : 2 * block]
    ).reshape(ngaussian, nstate)
    q = parameters[2 * block : 2 * block + ngaussian]
    p = parameters[2 * block + ngaussian :]
    return FrozenGaussianSpinorStateV251(
        q=q,
        p=p,
        widths=widths,
        coefficients=coefficients,
        time_au=float(time_au),
    ).validate(require_normalized=False)


def build_frozen_gaussian_spinor_matrices_v251(state, model):
    """Return exact overlap and Hamiltonian matrices for the quadratic contract."""

    state = state.validate(require_normalized=False)
    model = model.validate()
    if state.nstate != model.nstate:
        raise ValueError("state and quadratic Hamiltonian electronic dimensions differ.")
    dimension = state.ngaussian * state.nstate
    overlap = np.zeros((dimension, dimension), dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    identity = np.eye(state.nstate, dtype=complex)
    for i in range(state.ngaussian):
        si = slice(i * state.nstate, (i + 1) * state.nstate)
        for j in range(state.ngaussian):
            sj = slice(j * state.nstate, (j + 1) * state.nstate)
            moments = _cross_moments_v251(
                state.q[i],
                state.p[i],
                state.widths[i],
                state.q[j],
                state.p[j],
                state.widths[j],
                maximum_order=2,
            )
            kinetic = _kinetic_polynomial_v251(
                state.q[j],
                state.p[j],
                state.widths[j],
                model.mass_au,
            )
            overlap[si, sj] = moments[0] * identity
            hamiltonian[si, sj] = (
                _integrate_polynomial_v251(kinetic, moments) * identity
                + moments[0] * model.H0
                + moments[1] * model.H1
                + moments[2] * model.H2
            )
    if _hermiticity_residual_v251(overlap) > 2.0e-10:
        raise ValueError("analytic spinor overlap matrix is not Hermitian.")
    if _hermiticity_residual_v251(hamiltonian) > 2.0e-10:
        raise ValueError("analytic spinor Hamiltonian matrix is not Hermitian.")
    return overlap, hamiltonian


def variational_energy_v251(state, model):
    state = state.validate(require_normalized=False)
    overlap, hamiltonian = build_frozen_gaussian_spinor_matrices_v251(state, model)
    coefficients = state.coefficients.reshape(-1)
    norm = float(np.real(np.vdot(coefficients, overlap @ coefficients)))
    if norm <= 0.0:
        raise ValueError("cannot evaluate energy of a zero-norm variational state.")
    value = np.vdot(coefficients, hamiltonian @ coefficients) / norm
    if abs(float(np.imag(value))) > 2.0e-10:
        raise ValueError("variational energy has a non-negligible imaginary part.")
    return float(np.real(value))


@dataclass(frozen=True)
class MetricSolveReceiptV251:
    singular_values: np.ndarray
    cutoff: float
    rank: int
    nullity: int
    retained_condition_number: float
    rhs_norm: float
    null_rhs_norm: float
    null_rhs_relative: float
    linear_residual_norm: float
    linear_residual_relative: float
    velocity_norm: float

    def __post_init__(self):
        object.__setattr__(
            self,
            "singular_values",
            np.asarray(self.singular_values, dtype=float).copy(),
        )

    def validate(self, parameter_count=None):
        if self.singular_values.ndim != 1 or len(self.singular_values) < 1:
            raise ValueError("metric receipt requires a nonempty singular spectrum.")
        if not np.all(np.isfinite(self.singular_values)) or np.any(
            self.singular_values < 0.0
        ):
            raise ValueError("metric singular values must be finite and nonnegative.")
        if parameter_count is not None and len(self.singular_values) != int(
            parameter_count
        ):
            raise ValueError("metric spectrum length disagrees with parameter count.")
        if int(self.rank) != self.rank or int(self.nullity) != self.nullity:
            raise ValueError("metric rank and nullity must be integers.")
        if int(self.rank) < 1 or int(self.nullity) < 0:
            raise ValueError("metric rank/nullity are invalid.")
        if int(self.rank) + int(self.nullity) != len(self.singular_values):
            raise ValueError("metric rank and nullity do not span parameter space.")
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
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"metric receipt {name} must be finite and nonnegative.")
        if float(self.retained_condition_number) < 1.0:
            raise ValueError("retained metric condition number cannot be below one.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "singular_values": self.singular_values.tolist(),
            "cutoff": float(self.cutoff),
            "rank": int(self.rank),
            "nullity": int(self.nullity),
            "retained_condition_number": float(self.retained_condition_number),
            "rhs_norm": float(self.rhs_norm),
            "null_rhs_norm": float(self.null_rhs_norm),
            "null_rhs_relative": float(self.null_rhs_relative),
            "linear_residual_norm": float(self.linear_residual_norm),
            "linear_residual_relative": float(self.linear_residual_relative),
            "velocity_norm": float(self.velocity_norm),
        }


def solve_variational_metric_v251(
    metric,
    rhs,
    *,
    settings=VariationalMetricSettingsV251(),
):
    """Solve a positive-semidefinite McLachlan metric with compatible SVD nulls."""

    settings = settings.validate()
    metric = np.asarray(metric, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    if metric.ndim != 2 or metric.shape[0] < 1 or metric.shape[0] != metric.shape[1]:
        raise ValueError("variational metric must be a nonempty square matrix.")
    if rhs.shape != (metric.shape[0],):
        raise ValueError("variational metric RHS has incompatible shape.")
    if not np.all(np.isfinite(metric)) or not np.all(np.isfinite(rhs)):
        raise ValueError("variational metric system contains non-finite data.")
    symmetry = _scaled_norm_v251(metric, metric.T)
    if symmetry > settings.structural_tolerance:
        raise ValueError("variational metric must be real symmetric.")
    metric = 0.5 * (metric + metric.T)
    eigen_minimum = float(np.min(np.linalg.eigvalsh(metric)))
    scale = max(float(np.linalg.norm(metric, ord=2)), 1.0)
    if eigen_minimum < -settings.structural_tolerance * scale:
        raise ValueError("variational metric is not positive semidefinite.")

    left, singular_values, right_h = np.linalg.svd(metric, full_matrices=True)
    maximum = float(singular_values[0])
    cutoff = max(
        float(settings.metric_absolute_cutoff),
        float(settings.metric_relative_cutoff) * maximum,
    )
    retained = singular_values > cutoff
    rank = int(np.count_nonzero(retained))
    nullity = int(len(singular_values) - rank)
    if rank < 1:
        raise ValueError("variational metric has no retained direction.")
    minimum_retained = float(singular_values[retained][-1])
    condition = float(maximum / minimum_retained)
    if condition > float(settings.maximum_retained_condition_number):
        raise ValueError("retained variational metric is too ill-conditioned.")

    projected_rhs = left.T @ rhs
    null_rhs_norm = float(np.linalg.norm(projected_rhs[~retained]))
    rhs_norm = float(np.linalg.norm(rhs))
    rhs_scale = max(rhs_norm, 1.0e-30)
    null_rhs_relative = null_rhs_norm / rhs_scale
    if null_rhs_relative > float(settings.null_rhs_relative_tolerance):
        raise ValueError("variational metric RHS is incompatible with its null space.")

    inverse = np.zeros_like(singular_values)
    inverse[retained] = 1.0 / singular_values[retained]
    velocity = right_h.T @ (inverse * projected_rhs)
    residual_norm = float(np.linalg.norm(metric @ velocity - rhs))
    residual_relative = residual_norm / rhs_scale
    if residual_relative > float(settings.linear_residual_relative_tolerance):
        raise ValueError("SVD variational metric solve has an excessive residual.")
    receipt = MetricSolveReceiptV251(
        singular_values=singular_values,
        cutoff=cutoff,
        rank=rank,
        nullity=nullity,
        retained_condition_number=condition,
        rhs_norm=rhs_norm,
        null_rhs_norm=null_rhs_norm,
        null_rhs_relative=null_rhs_relative,
        linear_residual_norm=residual_norm,
        linear_residual_relative=residual_relative,
        velocity_norm=float(np.linalg.norm(velocity)),
    ).validate(parameter_count=len(rhs))
    return velocity, receipt


@dataclass(frozen=True)
class VariationalMetricSystemV251:
    metric: np.ndarray
    rhs: np.ndarray
    velocity: np.ndarray
    solve_receipt: MetricSolveReceiptV251
    generalized_norm: float
    energy_hartree: float
    settings: VariationalMetricSettingsV251

    def __post_init__(self):
        for name in ("metric", "rhs", "velocity"):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=float).copy()
            )

    def validate(self):
        settings = self.settings.validate()
        if self.metric.ndim != 2 or self.metric.shape[0] != self.metric.shape[1]:
            raise ValueError("stored variational metric must be square.")
        count = self.metric.shape[0]
        if self.rhs.shape != (count,) or self.velocity.shape != (count,):
            raise ValueError("stored variational metric vectors have incompatible shape.")
        if not all(
            np.all(np.isfinite(value))
            for value in (self.metric, self.rhs, self.velocity)
        ):
            raise ValueError("stored variational metric system is non-finite.")
        expected_velocity, expected_receipt = solve_variational_metric_v251(
            self.metric, self.rhs, settings=settings
        )
        if _scaled_norm_v251(self.velocity, expected_velocity) > settings.structural_tolerance:
            raise ValueError("stored TDVP velocity disagrees with the SVD metric solve.")
        if _scaled_norm_v251(
            self.solve_receipt.singular_values,
            expected_receipt.singular_values,
        ) > settings.structural_tolerance:
            raise ValueError("stored metric spectrum disagrees with the metric.")
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
            if abs(float(getattr(self.solve_receipt, name)) - float(getattr(expected_receipt, name))) > settings.structural_tolerance:
                raise ValueError(f"stored metric receipt {name} is inconsistent.")
        if self.solve_receipt.rank != expected_receipt.rank or self.solve_receipt.nullity != expected_receipt.nullity:
            raise ValueError("stored metric rank/nullity are inconsistent.")
        for name in ("generalized_norm", "energy_hartree"):
            if not np.isfinite(float(getattr(self, name))):
                raise ValueError(f"stored {name} must be finite.")
        if float(self.generalized_norm) <= 0.0:
            raise ValueError("stored generalized norm must be positive.")
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


def _tangent_terms_v251(state):
    """Return (Gaussian index, electronic vector, polynomial) in pack order."""

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
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray(
                    [
                        -state.widths[i] * state.q[i] - 1.0j * state.p[i],
                        state.widths[i],
                    ],
                    dtype=complex,
                ),
            )
        )
    for i in range(state.ngaussian):
        terms.append(
            (
                i,
                state.coefficients[i],
                np.asarray([-1.0j * state.q[i], 1.0j], dtype=complex),
            )
        )
    if len(terms) != state.parameter_count:
        raise AssertionError("v0.25.1 tangent count disagrees with parameter layout.")
    return terms


def build_variational_metric_system_v251(
    state,
    model,
    *,
    settings=VariationalMetricSettingsV251(),
):
    """Build and solve the exact fixed-width McLachlan metric system."""

    settings = settings.validate()
    state = state.validate(require_normalized=False)
    model = model.validate()
    if state.nstate != model.nstate:
        raise ValueError("state and model electronic dimensions differ.")
    terms = _tangent_terms_v251(state)
    count = len(terms)
    metric = np.zeros((count, count), dtype=float)
    rhs = np.zeros(count, dtype=float)
    moment_cache = {}
    for i in range(state.ngaussian):
        for j in range(state.ngaussian):
            moment_cache[(i, j)] = _cross_moments_v251(
                state.q[i],
                state.p[i],
                state.widths[i],
                state.q[j],
                state.p[j],
                state.widths[j],
                maximum_order=3,
            )

    for mu, (i, vector_i, polynomial_i) in enumerate(terms):
        conjugate_i = np.conj(polynomial_i)
        for nu, (j, vector_j, polynomial_j) in enumerate(terms):
            polynomial = np.convolve(conjugate_i, polynomial_j)
            nuclear = _integrate_polynomial_v251(
                polynomial, moment_cache[(i, j)]
            )
            metric[mu, nu] = float(np.real(np.vdot(vector_i, vector_j) * nuclear))

        projection = 0.0 + 0.0j
        for j in range(state.ngaussian):
            kinetic = _kinetic_polynomial_v251(
                state.q[j],
                state.p[j],
                state.widths[j],
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
    return VariationalMetricSystemV251(
        metric=metric,
        rhs=rhs,
        velocity=velocity,
        solve_receipt=solve_receipt,
        generalized_norm=state.generalized_norm,
        energy_hartree=variational_energy_v251(state, model),
        settings=settings,
    ).validate()


@dataclass(frozen=True)
class ImplicitMidpointTDVPStepV251:
    start: FrozenGaussianSpinorStateV251
    end: FrozenGaussianSpinorStateV251
    dt_au: float
    model: QuadraticSpinHamiltonianV251
    settings: VariationalMetricSettingsV251
    midpoint_parameters: np.ndarray
    midpoint_system: VariationalMetricSystemV251
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
    def energy_change_hartree(self):
        return float(self.end_energy_hartree - self.start_energy_hartree)

    @property
    def norm_change(self):
        return float(self.end_norm - self.start_norm)

    def validate(self):
        settings = self.settings.validate()
        model = self.model.validate()
        self.start.validate(
            require_normalized=True,
            tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
        )
        self.end.validate(
            require_normalized=False,
            tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
        )
        if self.start.nstate != model.nstate or self.end.nstate != model.nstate:
            raise ValueError("step state/model electronic dimensions differ.")
        if self.start.ngaussian != self.end.ngaussian:
            raise ValueError("Gaussian count changed inside a v0.25.1 step.")
        if not np.array_equal(self.start.widths, self.end.widths):
            raise ValueError("frozen Gaussian widths changed inside a v0.25.1 step.")
        self.end.validate(
            require_normalized=True,
            tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
        )
        dt = float(self.dt_au)
        if not np.isfinite(dt) or dt == 0.0:
            raise ValueError("implicit TDVP step dt must be finite and nonzero.")
        if abs(self.end.time_au - (self.start.time_au + dt)) > settings.structural_tolerance:
            raise ValueError("implicit TDVP endpoint time disagrees with dt.")
        theta_start = pack_variational_parameters_v251(self.start)
        theta_end = pack_variational_parameters_v251(self.end)
        expected_midpoint = 0.5 * (theta_start + theta_end)
        if _scaled_norm_v251(self.midpoint_parameters, expected_midpoint) > settings.structural_tolerance:
            raise ValueError("stored implicit midpoint parameters are inconsistent.")
        midpoint_state = state_from_variational_parameters_v251(
            expected_midpoint,
            widths=self.start.widths,
            nstate=self.start.nstate,
            time_au=self.start.time_au + 0.5 * dt,
        )
        expected_system = build_variational_metric_system_v251(
            midpoint_state, model, settings=settings
        )
        self.midpoint_system.validate()
        if self.midpoint_system.settings.as_dict() != settings.as_dict():
            raise ValueError("stored midpoint metric settings are inconsistent.")
        for name in ("metric", "rhs", "velocity"):
            if _scaled_norm_v251(
                getattr(self.midpoint_system, name), getattr(expected_system, name)
            ) > settings.structural_tolerance:
                raise ValueError(f"stored midpoint {name} is inconsistent.")
        if _scaled_norm_v251(
            self.midpoint_system.solve_receipt.singular_values,
            expected_system.solve_receipt.singular_values,
        ) > settings.structural_tolerance:
            raise ValueError("stored midpoint metric spectrum is inconsistent.")
        expected_residual = theta_end - theta_start - dt * expected_system.velocity
        if self.nonlinear_residual.shape != expected_residual.shape or _scaled_norm_v251(
            self.nonlinear_residual, expected_residual
        ) > settings.structural_tolerance:
            raise ValueError("stored nonlinear midpoint residual is inconsistent.")
        expected_residual_norm = float(np.linalg.norm(expected_residual))
        if abs(float(self.nonlinear_residual_norm) - expected_residual_norm) > settings.structural_tolerance:
            raise ValueError("stored nonlinear residual norm is inconsistent.")
        if expected_residual_norm > settings.nonlinear_residual_tolerance:
            raise ValueError("implicit midpoint nonlinear residual exceeds tolerance.")
        if type(self.nonlinear_success) is not bool or not self.nonlinear_success:
            raise ValueError("implicit midpoint nonlinear solver did not report success.")
        if not isinstance(self.nonlinear_status, (int, np.integer)) or int(self.nonlinear_status) <= 0:
            raise ValueError("implicit midpoint nonlinear status is not successful.")
        if not isinstance(self.nonlinear_message, str) or not self.nonlinear_message:
            raise ValueError("implicit midpoint nonlinear message is missing.")
        if (
            isinstance(self.nonlinear_function_evaluations, (bool, np.bool_))
            or not isinstance(
                self.nonlinear_function_evaluations, (int, np.integer)
            )
            or int(self.nonlinear_function_evaluations) < 1
        ):
            raise ValueError("implicit midpoint function-evaluation count is invalid.")
        if self.nonlinear_jacobian_evaluations is not None:
            if (
                isinstance(self.nonlinear_jacobian_evaluations, (bool, np.bool_))
                or not isinstance(
                    self.nonlinear_jacobian_evaluations, (int, np.integer)
                )
                or int(self.nonlinear_jacobian_evaluations) < 0
            ):
                raise ValueError(
                    "implicit midpoint Jacobian-evaluation count is invalid."
                )
        if not np.isfinite(float(self.predictor_residual_norm)) or float(
            self.predictor_residual_norm
        ) < 0.0:
            raise ValueError("implicit midpoint predictor residual is invalid.")
        expected_start_norm = self.start.generalized_norm
        expected_end_norm = self.end.generalized_norm
        expected_start_energy = variational_energy_v251(self.start, model)
        expected_end_energy = variational_energy_v251(self.end, model)
        for name, stored, expected in (
            ("start_norm", self.start_norm, expected_start_norm),
            ("end_norm", self.end_norm, expected_end_norm),
            ("start_energy_hartree", self.start_energy_hartree, expected_start_energy),
            ("end_energy_hartree", self.end_energy_hartree, expected_end_energy),
        ):
            if abs(float(stored) - float(expected)) > settings.structural_tolerance:
                raise ValueError(f"stored {name} is inconsistent.")
        if abs(self.norm_change) > settings.maximum_step_norm_drift:
            raise ValueError("implicit TDVP step norm drift exceeds the release gate.")
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
        }


def implicit_midpoint_tdvp_step_v251(
    state,
    model,
    dt_au,
    *,
    settings=VariationalMetricSettingsV251(),
):
    """Advance one signed, fully implicit frozen-width TDVP step."""

    settings = settings.validate()
    model = model.validate()
    state = state.validate(
        require_normalized=True,
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
    )
    if state.nstate != model.nstate:
        raise ValueError("state and model electronic dimensions differ.")
    dt = float(dt_au)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError("implicit TDVP step dt must be finite and nonzero.")
    theta_start = pack_variational_parameters_v251(state)
    initial_system = build_variational_metric_system_v251(
        state, model, settings=settings
    )
    predictor = theta_start + dt * initial_system.velocity

    def residual(theta_end):
        midpoint = 0.5 * (theta_start + np.asarray(theta_end, dtype=float))
        midpoint_state = state_from_variational_parameters_v251(
            midpoint,
            widths=state.widths,
            nstate=state.nstate,
            time_au=state.time_au + 0.5 * dt,
        )
        midpoint_system = build_variational_metric_system_v251(
            midpoint_state, model, settings=settings
        )
        return np.asarray(theta_end, dtype=float) - theta_start - dt * midpoint_system.velocity

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
            "implicit midpoint TDVP solve failed: "
            f"success={bool(solution.success)}, status={int(solution.status)}, "
            f"residual={final_residual_norm:.6e}, message={solution.message}"
        )
    end = state_from_variational_parameters_v251(
        theta_end,
        widths=state.widths,
        nstate=state.nstate,
        time_au=state.time_au + dt,
    ).validate(
        require_normalized=True,
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
    )
    midpoint_parameters = 0.5 * (theta_start + theta_end)
    midpoint_state = state_from_variational_parameters_v251(
        midpoint_parameters,
        widths=state.widths,
        nstate=state.nstate,
        time_au=state.time_au + 0.5 * dt,
    )
    midpoint_system = build_variational_metric_system_v251(
        midpoint_state, model, settings=settings
    )
    return ImplicitMidpointTDVPStepV251(
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
        start_energy_hartree=variational_energy_v251(state, model),
        end_energy_hartree=variational_energy_v251(end, model),
    ).validate()


V251_TDVP_CLAIMS = {
    "frozen_width_multigaussian_tdvp_metric_validated": True,
    "implicit_midpoint_nonlinear_solve_validated": True,
    "svd_metric_rank_and_compatible_nullspace_validated": True,
    "complete_spinor_quadratic_soc_validated": True,
    "single_packet_harmonic_continuous_reduction_validated": True,
    "gaussian_permutation_covariance_validated": True,
    "constant_electronic_gauge_covariance_validated": True,
    "adaptive_gaussian_width_tdvp_validated": False,
    "dynamic_spawning_validated": False,
    "dynamic_pruning_validated": False,
    "coordinate_dependent_electronic_gauge_covariance_validated": False,
    "multidimensional_multigaussian_tdvp_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
    "general_ab_initio_soc_dynamics_accuracy_validated": False,
}


@dataclass(frozen=True)
class FrozenWidthMultiGaussianTrajectoryV251:
    initial_state: FrozenGaussianSpinorStateV251
    final_state: FrozenGaussianSpinorStateV251
    model: QuadraticSpinHamiltonianV251
    settings: VariationalMetricSettingsV251
    steps: tuple
    claims: dict = field(default_factory=lambda: dict(V251_TDVP_CLAIMS))

    @property
    def maximum_norm_drift(self):
        initial = self.initial_state.generalized_norm
        values = [abs(step.end_norm - initial) for step in self.steps]
        return float(max(values, default=0.0))

    @property
    def maximum_absolute_energy_drift_hartree(self):
        initial = variational_energy_v251(self.initial_state, self.model)
        values = [abs(step.end_energy_hartree - initial) for step in self.steps]
        return float(max(values, default=0.0))

    @property
    def minimum_metric_rank(self):
        values = [step.midpoint_system.solve_receipt.rank for step in self.steps]
        return int(min(values, default=self.initial_state.parameter_count))

    @property
    def maximum_metric_nullity(self):
        values = [step.midpoint_system.solve_receipt.nullity for step in self.steps]
        return int(max(values, default=0))

    @property
    def maximum_nonlinear_residual(self):
        return float(max((step.nonlinear_residual_norm for step in self.steps), default=0.0))

    def validate(self):
        model = self.model.validate()
        settings = self.settings.validate()
        self.initial_state.validate(
            require_normalized=True,
            tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
        )
        self.final_state.validate(
            require_normalized=True,
            tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift),
        )
        if self.initial_state.nstate != model.nstate:
            raise ValueError("trajectory state/model dimensions differ.")
        if type(self.claims) is not dict or any(
            type(value) is not bool for value in self.claims.values()
        ):
            raise TypeError("every v0.25.1 TDVP claim must be a native Boolean.")
        if self.claims != V251_TDVP_CLAIMS:
            raise ValueError("v0.25.1 TDVP claims differ from the frozen boundary.")
        previous = self.initial_state
        for step in self.steps:
            step.validate()
            if step.model.fingerprint() != model.fingerprint():
                raise ValueError("trajectory step model identity changed.")
            if step.settings.as_dict() != settings.as_dict():
                raise ValueError("trajectory step settings changed.")
            if _scaled_norm_v251(
                pack_variational_parameters_v251(step.start),
                pack_variational_parameters_v251(previous),
            ) > settings.structural_tolerance or abs(
                step.start.time_au - previous.time_au
            ) > settings.structural_tolerance:
                raise ValueError("v0.25.1 TDVP step chain is discontinuous.")
            previous = step.end
        if _scaled_norm_v251(
            pack_variational_parameters_v251(previous),
            pack_variational_parameters_v251(self.final_state),
        ) > settings.structural_tolerance or abs(
            previous.time_au - self.final_state.time_au
        ) > settings.structural_tolerance:
            raise ValueError("v0.25.1 TDVP final state is not bound to the last step.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIGAUSSIAN_TDVP_SCHEMA_V251,
            "ansatz": MULTIGAUSSIAN_TDVP_ANSATZ_V251,
            "variational_principle": VARIATIONAL_PRINCIPLE_V251,
            "integrator": VARIATIONAL_INTEGRATOR_V251,
            "metric_solver": VARIATIONAL_METRIC_SOLVER_V251,
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
            "claims": dict(self.claims),
        }

    def fingerprint(self):
        return _sha256_v251(self.as_dict())


def run_frozen_width_multigaussian_tdvp_v251(
    initial_state,
    model,
    *,
    dt_au,
    steps,
    settings=VariationalMetricSettingsV251(),
):
    settings = settings.validate()
    model = model.validate()
    dt = float(dt_au)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("forward multi-Gaussian TDVP dt must be positive.")
    if (
        isinstance(steps, (bool, np.bool_))
        or not isinstance(steps, (int, np.integer))
        or int(steps) < 0
    ):
        raise ValueError("multi-Gaussian TDVP steps must be a nonnegative integer.")
    state = initial_state.normalized()
    initial = state
    receipts = []
    for _ in range(int(steps)):
        receipt = implicit_midpoint_tdvp_step_v251(
            state, model, dt, settings=settings
        )
        receipts.append(receipt)
        state = receipt.end
    return FrozenWidthMultiGaussianTrajectoryV251(
        initial_state=initial,
        final_state=state,
        model=model,
        settings=settings,
        steps=tuple(receipts),
        claims=dict(V251_TDVP_CLAIMS),
    ).validate()


def reverse_frozen_width_multigaussian_tdvp_v251(trajectory):
    trajectory = trajectory.validate()
    state = trajectory.final_state
    receipts = []
    for forward_step in reversed(trajectory.steps):
        receipt = implicit_midpoint_tdvp_step_v251(
            state,
            trajectory.model,
            -forward_step.dt_au,
            settings=trajectory.settings,
        )
        receipts.append(receipt)
        state = receipt.end
    return FrozenWidthMultiGaussianTrajectoryV251(
        initial_state=trajectory.final_state,
        final_state=state,
        model=trajectory.model,
        settings=trajectory.settings,
        steps=tuple(receipts),
        claims=dict(V251_TDVP_CLAIMS),
    ).validate()
