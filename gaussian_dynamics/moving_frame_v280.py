"""Flat coordinate-dependent electronic frames for Gaussian dynamics v0.28.0.

The admitted v0.28.0 manifold is

    Psi(R) = sum_I g_I(R) Phi(R) W(R,q_I) c_I,

where ``c_I`` is stored in the electronic frame at the packet centre and
``W(R,q)=G(R)^† G(q)`` is exact parallel transport for a flat pure-gauge frame
``Phi(R)=Phi_ref G(R)``. The physical fixed-reference coefficient is therefore
``G(q_I)c_I``.

This first v0.28 milestone deliberately admits only analytically trivializable flat
connections. Nonzero curvature and live molecular-SOC trajectories remain closed.
"""

from dataclasses import dataclass
import numpy as np

from .correlated_gaussian_tdvp_v270 import (
    CorrelatedGaussianSpinorStateV270,
    CorrelatedVariationalSettingsV270,
    _pack_velocity_blocks_v270,
    _velocity_blocks_v270,
    build_correlated_metric_system_v270,
    correlated_implicit_midpoint_step_v270,
    evaluate_correlated_state_v270,
)

MOVING_FRAME_SCHEMA_V280 = "gnd-flat-moving-electronic-frame-v0.28.0"
MOVING_STATE_SCHEMA_V280 = "gnd-moving-frame-correlated-state-v0.28.0"
CONNECTION_CONVENTION_V280 = "D_a=G^dagger partial_a G; covariant derivative partial_a+D_a"
TRANSPORT_CONVENTION_V280 = "W(R,q)=G(R)^dagger G(q)"
CLAIM_BOUNDARY_V280 = (
    "analytic flat pure-gauge electronic frames exactly trivializable to the sealed "
    "v0.27.0 fixed-frame correlated-Gaussian manifold"
)


def _scaled_norm_v280(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right) / scale)


def _apply_matrix_v280(matrix, vector):
    return np.einsum("...ab,...b->...a", matrix, vector, optimize=True)


@dataclass(frozen=True)
class FlatMovingFrameV280:
    generator: np.ndarray
    phase_gradient: np.ndarray
    phase_hessian: np.ndarray
    phase_offset: float = 0.0
    right_unitary: np.ndarray | None = None
    label: str = "analytic flat moving frame"

    def __post_init__(self):
        generator = np.asarray(self.generator, dtype=complex).copy()
        gradient = np.asarray(self.phase_gradient, dtype=float).copy()
        hessian = np.asarray(self.phase_hessian, dtype=float).copy()
        if self.right_unitary is None:
            right = np.eye(generator.shape[0], dtype=complex) if generator.ndim == 2 else np.asarray([], dtype=complex)
        else:
            right = np.asarray(self.right_unitary, dtype=complex).copy()
        object.__setattr__(self, "generator", generator)
        object.__setattr__(self, "phase_gradient", gradient)
        object.__setattr__(self, "phase_hessian", hessian)
        object.__setattr__(self, "right_unitary", right)

    @property
    def nstate(self):
        return int(self.generator.shape[0]) if self.generator.ndim == 2 else 0

    @property
    def ndim(self):
        return int(self.phase_gradient.shape[0]) if self.phase_gradient.ndim == 1 else 0

    def validate(self, tolerance=2.0e-11):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("moving-frame tolerance must be finite and positive.")
        if self.generator.ndim != 2 or self.nstate < 1 or self.generator.shape != (self.nstate, self.nstate):
            raise ValueError("moving-frame generator must be a nonempty square matrix.")
        if not np.all(np.isfinite(self.generator)) or _scaled_norm_v280(self.generator, self.generator.conj().T) > tolerance:
            raise ValueError("moving-frame generator must be finite and Hermitian.")
        if self.phase_gradient.ndim != 1 or self.ndim < 1:
            raise ValueError("phase gradient must be a nonempty vector.")
        if self.phase_hessian.shape != (self.ndim, self.ndim):
            raise ValueError("phase Hessian has an incompatible shape.")
        if not np.all(np.isfinite(self.phase_gradient)) or not np.all(np.isfinite(self.phase_hessian)):
            raise ValueError("moving-frame phase data must be finite.")
        if _scaled_norm_v280(self.phase_hessian, self.phase_hessian.T) > tolerance:
            raise ValueError("phase Hessian must be symmetric.")
        if not np.isfinite(float(self.phase_offset)):
            raise ValueError("phase offset must be finite.")
        if self.right_unitary.shape != (self.nstate, self.nstate):
            raise ValueError("right-unitary gauge has an incompatible shape.")
        if _scaled_norm_v280(self.right_unitary.conj().T @ self.right_unitary, np.eye(self.nstate)) > tolerance:
            raise ValueError("right-unitary gauge must be unitary.")
        if not str(self.label):
            raise ValueError("moving-frame label cannot be empty.")
        return self

    def theta(self, coordinates):
        self.validate()
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape == () or coordinates.shape[-1] != self.ndim or not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates must be finite with final axis ndim.")
        return (
            float(self.phase_offset)
            + np.einsum("...a,a->...", coordinates, self.phase_gradient, optimize=True)
            + 0.5 * np.einsum("...a,ab,...b->...", coordinates, self.phase_hessian, coordinates, optimize=True)
        )

    def theta_gradient(self, coordinates):
        self.validate()
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape == () or coordinates.shape[-1] != self.ndim or not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates must be finite with final axis ndim.")
        return self.phase_gradient + np.einsum("ab,...b->...a", self.phase_hessian, coordinates, optimize=True)

    @property
    def effective_generator(self):
        self.validate()
        return self.right_unitary.conj().T @ self.generator @ self.right_unitary

    def unitary(self, coordinates):
        self.validate()
        theta = np.asarray(self.theta(coordinates), dtype=float)
        values, vectors = np.linalg.eigh(self.generator)
        phases = np.exp(1.0j * theta[..., None] * values)
        exponential = np.einsum("ai,...i,bi->...ab", vectors, phases, vectors.conj(), optimize=True)
        return exponential @ self.right_unitary

    def connection(self, coordinates):
        self.validate()
        gradient = self.theta_gradient(coordinates)
        return 1.0j * gradient[..., :, None, None] * self.effective_generator

    def curvature(self, coordinates):
        coordinates = np.asarray(coordinates, dtype=float)
        self.theta(coordinates)
        shape = coordinates.shape[:-1] + (self.ndim, self.ndim, self.nstate, self.nstate)
        return np.zeros(shape, dtype=complex)

    def transporter(self, coordinates, center):
        self.validate()
        coordinates = np.asarray(coordinates, dtype=float)
        center = np.asarray(center, dtype=float)
        if center.shape != (self.ndim,) or not np.all(np.isfinite(center)):
            raise ValueError("transport center must be a finite ndim vector.")
        g_r = self.unitary(coordinates)
        g_q = self.unitary(center)
        return np.einsum("...ba,bc->...ac", g_r.conj(), g_q, optimize=True)

    def gauge_rotated(self, unitary):
        unitary = np.asarray(unitary, dtype=complex)
        if unitary.shape != (self.nstate, self.nstate) or _scaled_norm_v280(unitary.conj().T @ unitary, np.eye(self.nstate)) > 2.0e-11:
            raise ValueError("moving-frame gauge rotation must be unitary.")
        return FlatMovingFrameV280(
            self.generator,
            self.phase_gradient,
            self.phase_hessian,
            self.phase_offset,
            self.right_unitary @ unitary,
            self.label + " [gauge rotated]",
        ).validate()


def require_flat_moving_frame_v280(frame, *, curvature=None, trivialization_available=True, tolerance=2.0e-10):
    """Fail closed unless the connection is flat and exactly trivializable."""
    if type(trivialization_available) is not bool:
        raise TypeError("trivialization availability must be a native Boolean.")
    if not trivialization_available:
        raise ValueError("v0.28.0 requires an exact electronic-frame trivialization.")
    frame = frame.validate()
    tolerance = float(tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("flatness tolerance must be finite and positive.")
    if curvature is None:
        curvature = frame.curvature(np.zeros(frame.ndim, dtype=float))
    curvature = np.asarray(curvature, dtype=complex)
    expected = (frame.ndim, frame.ndim, frame.nstate, frame.nstate)
    if curvature.shape != expected or not np.all(np.isfinite(curvature)):
        raise ValueError("curvature tensor has an invalid shape or values.")
    if float(np.linalg.norm(curvature)) > tolerance:
        raise ValueError("v0.28.0 rejects non-flat electronic connections.")
    return frame


@dataclass(frozen=True)
class MovingFrameCorrelatedStateV280:
    q: np.ndarray
    p: np.ndarray
    width_matrices: np.ndarray
    chirp_matrices: np.ndarray
    center_coefficients: np.ndarray
    time_au: float = 0.0

    def __post_init__(self):
        for name in ("q", "p", "width_matrices", "chirp_matrices"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())
        object.__setattr__(self, "center_coefficients", np.asarray(self.center_coefficients, dtype=complex).copy())

    @property
    def ngaussian(self):
        return int(self.q.shape[0]) if self.q.ndim == 2 else 0

    @property
    def ndim(self):
        return int(self.q.shape[1]) if self.q.ndim == 2 else 0

    @property
    def nstate(self):
        return int(self.center_coefficients.shape[1]) if self.center_coefficients.ndim == 2 else 0

    def validate(self, frame, *, require_normalized=False, tolerance=5.0e-10):
        frame = frame.validate()
        if self.q.ndim != 2 or self.ngaussian < 1 or self.ndim != frame.ndim:
            raise ValueError("moving state q has an invalid shape or frame dimension.")
        if self.p.shape != self.q.shape:
            raise ValueError("moving state q and p must have identical shapes.")
        if self.width_matrices.shape != (self.ngaussian, self.ndim, self.ndim):
            raise ValueError("moving-state widths have an invalid shape.")
        if self.chirp_matrices.shape != self.width_matrices.shape:
            raise ValueError("moving-state chirps have an invalid shape.")
        if self.center_coefficients.shape != (self.ngaussian, frame.nstate):
            raise ValueError("moving-state coefficients have an invalid shape.")
        moving_to_fixed_state_v280(self, frame).validate(require_normalized=require_normalized, tolerance=tolerance)
        return self

    def gauge_rotated(self, unitary):
        unitary = np.asarray(unitary, dtype=complex)
        if unitary.shape != (self.nstate, self.nstate) or _scaled_norm_v280(unitary.conj().T @ unitary, np.eye(self.nstate)) > 2.0e-11:
            raise ValueError("state gauge rotation must be unitary.")
        coefficients = np.einsum("ab,Ib->Ia", unitary.conj().T, self.center_coefficients, optimize=True)
        return MovingFrameCorrelatedStateV280(
            self.q, self.p, self.width_matrices, self.chirp_matrices, coefficients, self.time_au
        )


def moving_to_fixed_state_v280(state, frame):
    frame = frame.validate()
    if not isinstance(state, MovingFrameCorrelatedStateV280):
        raise TypeError("moving-to-fixed conversion requires a v0.28 moving state.")
    if state.q.ndim != 2 or state.ndim != frame.ndim or state.nstate != frame.nstate:
        raise ValueError("moving state and frame dimensions are incompatible.")
    coefficients = _apply_matrix_v280(frame.unitary(state.q), state.center_coefficients)
    return CorrelatedGaussianSpinorStateV270(
        state.q, state.p, state.width_matrices, state.chirp_matrices, coefficients, state.time_au
    )


def fixed_to_moving_state_v280(state, frame):
    frame = frame.validate()
    state = state.validate(require_normalized=False)
    if state.ndim != frame.ndim or state.nstate != frame.nstate:
        raise ValueError("fixed state and moving frame dimensions are incompatible.")
    g_q = frame.unitary(state.q)
    center = np.einsum("...ba,...b->...a", g_q.conj(), state.coefficients, optimize=True)
    return MovingFrameCorrelatedStateV280(
        state.q, state.p, state.width_matrices, state.chirp_matrices, center, state.time_au
    )


def moving_frame_hamiltonian_v280(model, frame, coordinates):
    model = model.validate()
    frame = frame.validate()
    if model.ndim != frame.ndim or model.nstate != frame.nstate:
        raise ValueError("model and moving frame dimensions are incompatible.")
    h_ref = model.hamiltonian(coordinates)
    g = frame.unitary(coordinates)
    return np.einsum("...ba,...bc,...cd->...ad", g.conj(), h_ref, g, optimize=True)


def evaluate_moving_section_v280(state, frame, points):
    state = state.validate(frame, require_normalized=False)
    points = np.asarray(points, dtype=float)
    if points.shape == () or points.shape[-1] != state.ndim or not np.all(np.isfinite(points)):
        raise ValueError("evaluation points must be finite with final axis ndim.")
    section = np.zeros(points.shape[:-1] + (state.nstate,), dtype=complex)
    for packet in range(state.ngaussian):
        displacement = points - state.q[packet]
        width = state.width_matrices[packet]
        chirp = state.chirp_matrices[packet]
        logdet = float(np.sum(np.log(np.linalg.eigvalsh(width))))
        normalization = float(np.exp(0.25 * (logdet - state.ndim * np.log(np.pi))))
        exponent = (
            -0.5 * np.einsum("...a,ab,...b->...", displacement, width, displacement)
            + 0.5j * np.einsum("...a,ab,...b->...", displacement, chirp, displacement)
            + 1.0j * np.einsum("...a,a->...", displacement, state.p[packet])
        )
        gaussian = normalization * np.exp(exponent)
        transported = _apply_matrix_v280(frame.transporter(points, state.q[packet]), state.center_coefficients[packet])
        section += gaussian[..., None] * transported
    return section


def evaluate_moving_physical_v280(state, frame, points):
    return _apply_matrix_v280(frame.unitary(points), evaluate_moving_section_v280(state, frame, points))


def moving_frame_velocity_v280(state, model, frame, *, settings=CorrelatedVariationalSettingsV270(), active_shape_mask=None):
    state = state.validate(frame, require_normalized=False)
    fixed = moving_to_fixed_state_v280(state, frame).validate(require_normalized=False)
    system = build_correlated_metric_system_v270(fixed, model, settings=settings, active_shape_mask=active_shape_mask)
    cfdot, qdot, pdot, edot, bdot = _velocity_blocks_v270(fixed, system.velocity)
    g_q = frame.unitary(state.q)
    connection = frame.connection(state.q)
    cmdot = np.empty_like(state.center_coefficients)
    for packet in range(state.ngaussian):
        local_fixed_dot = g_q[packet].conj().T @ cfdot[packet]
        contracted = np.einsum("a,aij->ij", qdot[packet], connection[packet], optimize=True)
        cmdot[packet] = local_fixed_dot - contracted @ state.center_coefficients[packet]
    return _pack_velocity_blocks_v270(cmdot, qdot, pdot, edot, bdot), system


@dataclass(frozen=True)
class MovingFrameImplicitMidpointStepV280:
    start: MovingFrameCorrelatedStateV280
    end: MovingFrameCorrelatedStateV280
    fixed_step: object
    frame: FlatMovingFrameV280

    def validate(self):
        self.frame.validate()
        self.start.validate(self.frame, require_normalized=False)
        self.end.validate(self.frame, require_normalized=False)
        self.fixed_step.validate()
        fixed_start = moving_to_fixed_state_v280(self.start, self.frame)
        fixed_end = moving_to_fixed_state_v280(self.end, self.frame)
        if _scaled_norm_v280(fixed_start.q, self.fixed_step.start.q) > 1.0e-12:
            raise ValueError("moving-step start differs from the fixed trivialization.")
        if _scaled_norm_v280(fixed_start.coefficients, self.fixed_step.start.coefficients) > 2.0e-11:
            raise ValueError("moving-step start coefficients differ from the fixed trivialization.")
        if _scaled_norm_v280(fixed_end.coefficients, self.fixed_step.end.coefficients) > 2.0e-11:
            raise ValueError("moving-step endpoint differs from the fixed trivialization.")
        return self


def moving_frame_implicit_midpoint_step_v280(state, model, frame, dt_au, *, settings=CorrelatedVariationalSettingsV270(), active_shape_mask=None):
    state = state.validate(frame, require_normalized=False)
    fixed_start = moving_to_fixed_state_v280(state, frame).validate(require_normalized=False)
    fixed_step = correlated_implicit_midpoint_step_v270(
        fixed_start, model, dt_au, settings=settings, active_shape_mask=active_shape_mask
    )
    end = fixed_to_moving_state_v280(fixed_step.end, frame)
    return MovingFrameImplicitMidpointStepV280(state, end, fixed_step, frame).validate()


@dataclass(frozen=True)
class MovingFrameBasisEventV280:
    before: MovingFrameCorrelatedStateV280
    after: MovingFrameCorrelatedStateV280
    fixed_event: object
    frame: FlatMovingFrameV280

    @property
    def event_type(self):
        return str(self.fixed_event.event_kind)

    @property
    def reason(self):
        return str(self.fixed_event.reason)

    def validate(self):
        self.frame.validate()
        self.before.validate(self.frame, require_normalized=False)
        self.after.validate(self.frame, require_normalized=False)
        self.fixed_event.validate()
        fixed_before = moving_to_fixed_state_v280(self.before, self.frame)
        fixed_after = moving_to_fixed_state_v280(self.after, self.frame)
        if _scaled_norm_v280(fixed_before.coefficients, self.fixed_event.before.coefficients) > 2.0e-11:
            raise ValueError("moving basis-event input differs from its fixed trivialization.")
        if _scaled_norm_v280(fixed_after.coefficients, self.fixed_event.after.coefficients) > 2.0e-11:
            raise ValueError("moving basis-event output differs from its fixed trivialization.")
        return self


def adapt_moving_frame_basis_once_v280(state, model, frame, **kwargs):
    from .correlated_basis_adaptation_v270 import adapt_multidimensional_basis_once_v270
    state = state.validate(frame, require_normalized=False)
    fixed = moving_to_fixed_state_v280(state, frame)
    event = adapt_multidimensional_basis_once_v270(fixed, model, **kwargs)
    return MovingFrameBasisEventV280(
        fixed_to_moving_state_v280(event.before, frame),
        fixed_to_moving_state_v280(event.after, frame),
        event,
        frame,
    ).validate()


def reference_wavefunction_error_v280(state, frame, points):
    moving = evaluate_moving_physical_v280(state, frame, points)
    fixed = evaluate_correlated_state_v270(moving_to_fixed_state_v280(state, frame), points)
    return float(np.max(np.abs(moving - fixed), initial=0.0))


__all__ = [
    "CLAIM_BOUNDARY_V280", "CONNECTION_CONVENTION_V280", "TRANSPORT_CONVENTION_V280",
    "MOVING_FRAME_SCHEMA_V280", "MOVING_STATE_SCHEMA_V280", "FlatMovingFrameV280",
    "MovingFrameCorrelatedStateV280", "MovingFrameImplicitMidpointStepV280",
    "MovingFrameBasisEventV280", "adapt_moving_frame_basis_once_v280",
    "evaluate_moving_physical_v280", "evaluate_moving_section_v280",
    "fixed_to_moving_state_v280", "moving_frame_hamiltonian_v280",
    "moving_frame_implicit_midpoint_step_v280", "moving_frame_velocity_v280",
    "moving_to_fixed_state_v280", "reference_wavefunction_error_v280",
    "require_flat_moving_frame_v280",
]
