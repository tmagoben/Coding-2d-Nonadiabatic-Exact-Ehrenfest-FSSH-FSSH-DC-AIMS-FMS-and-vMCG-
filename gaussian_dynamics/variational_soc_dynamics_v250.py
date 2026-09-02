"""Symmetric time-dependent-variational SOC propagation for v0.25.0.

The validated ansatz in this release is the canonical single-nuclear-packet
Ehrenfest limit of a time-dependent variational principle,

    |Psi(t)> = |g(q(t),p(t))> sum_I c_I(t) |Phi_I(q(t))>.

It is not the full coupled multi-Gaussian TDVP.  Canonical nuclear variables use a
time-reversible kick--drift--kick update only when the generalized mass is constant.
Electronic amplitudes use endpoint Strang propagation and the unitary polar factor
of the cross-geometry overlap.  The polar factor is computed and certified by SVD;
"polar decomposition" and "SVD" are therefore complementary, not competing,
choices.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np

from .finite_manifold_transport_v233 import (
    FiniteManifoldOverlapPolicyV233,
    certified_transport_from_overlap_v233,
)
from .temporal_electronic import hermitian_exponential


VARIATIONAL_SOC_SCHEMA_V250 = "gnd-symmetric-variational-soc-trajectory-v0.25.0"
RESTRICTED_TDVP_ANSATZ_V250 = (
    "single canonical nuclear packet with a complete electronic spinor"
)
RESTRICTED_NUCLEAR_INTEGRATOR_V250 = (
    "symmetric velocity-Verlet for constant-mass canonical variables"
)
GENERAL_TDVP_INTEGRATOR_V250 = (
    "implicit midpoint/discrete variational solve for the coupled noncanonical manifold"
)
ELECTRONIC_INTEGRATOR_V250 = (
    "endpoint-Hamiltonian Strang step with right-to-left unitary polar transport"
)
POLAR_ALGORITHM_V250 = "SVD O=U Sigma V^dagger; W=U V^dagger"


def _canonical_v250(value):
    if isinstance(value, np.generic):
        return _canonical_v250(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_v250(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_v250(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v250(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("v0.25.0 canonical data cannot contain non-finite values.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.25.0 canonical value: {type(value).__name__}")


def _sha256_v250(value):
    payload = json.dumps(
        _canonical_v250(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_v250(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _scaled_frobenius_v250(left, right):
    left = np.asarray(left, dtype=complex)
    right = np.asarray(right, dtype=complex)
    if left.shape != right.shape:
        return float("inf")
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    return float(np.linalg.norm(left - right, ord="fro") / scale)


def _state_distance_v250(left, right):
    return max(
        float(np.max(np.abs(left.q - right.q))),
        float(np.max(np.abs(left.p - right.p))),
        float(np.max(np.abs(left.electronic_coefficients - right.electronic_coefficients))),
        abs(float(left.time_au) - float(right.time_au)),
    )


@dataclass(frozen=True)
class VariationalSOCIntegratorSettingsV250:
    """Frozen scope and tolerances for the v0.25.0 symmetric propagator."""

    overlap_policy: FiniteManifoldOverlapPolicyV233 = field(
        default_factory=lambda: FiniteManifoldOverlapPolicyV233(
            minimum_retained_singular_value=0.9,
            maximum_condition_number=10.0,
            maximum_principal_angle_radians=float(np.arccos(0.9)),
        )
    )
    mass_relative_tolerance: float = 1.0e-10
    mass_absolute_tolerance: float = 1.0e-12
    structural_tolerance: float = 1.0e-10
    step_binding_tolerance: float = 2.0e-9
    full_multi_gaussian_tdvp: bool = False
    adaptive_gaussian_widths: bool = False
    coordinate_dependent_mass: bool = False
    nuclear_integrator: str = RESTRICTED_NUCLEAR_INTEGRATOR_V250
    electronic_integrator: str = ELECTRONIC_INTEGRATOR_V250
    polar_algorithm: str = POLAR_ALGORITHM_V250

    def validate(self):
        self.overlap_policy.validate()
        for name in (
            "mass_relative_tolerance",
            "mass_absolute_tolerance",
            "structural_tolerance",
            "step_binding_tolerance",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        for name in (
            "full_multi_gaussian_tdvp",
            "adaptive_gaussian_widths",
            "coordinate_dependent_mass",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if self.full_multi_gaussian_tdvp:
            raise ValueError(
                "v0.25.0 does not admit full multi-Gaussian TDVP; use a future "
                "implicit-midpoint/discrete-variational implementation."
            )
        if self.adaptive_gaussian_widths:
            raise ValueError("v0.25.0 does not admit adaptive Gaussian-width dynamics.")
        if self.coordinate_dependent_mass:
            raise ValueError(
                "velocity-Verlet is not admitted for coordinate-dependent generalized mass."
            )
        if self.nuclear_integrator != RESTRICTED_NUCLEAR_INTEGRATOR_V250:
            raise ValueError("the v0.25.0 restricted nuclear integrator is frozen.")
        if self.electronic_integrator != ELECTRONIC_INTEGRATOR_V250:
            raise ValueError("the v0.25.0 electronic integrator is frozen.")
        if self.polar_algorithm != POLAR_ALGORITHM_V250:
            raise ValueError("the v0.25.0 SVD-polar algorithm is frozen.")
        return self

    def as_dict(self):
        payload = asdict(self)
        return _canonical_v250(payload)


@dataclass(frozen=True)
class CanonicalVariationalSOCStateV250:
    q: np.ndarray
    p: np.ndarray
    electronic_coefficients: np.ndarray
    time_au: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "q", np.asarray(self.q, dtype=float).copy())
        object.__setattr__(self, "p", np.asarray(self.p, dtype=float).copy())
        object.__setattr__(
            self,
            "electronic_coefficients",
            np.asarray(self.electronic_coefficients, dtype=complex).copy(),
        )

    @property
    def electronic_norm(self):
        return float(np.real(np.vdot(self.electronic_coefficients, self.electronic_coefficients)))

    def validate(self, tolerance=1.0e-10):
        if self.q.ndim != 1 or len(self.q) < 1:
            raise ValueError("variational SOC coordinates must be a nonempty vector.")
        if self.p.shape != self.q.shape:
            raise ValueError("variational SOC momentum and coordinate shapes differ.")
        if self.electronic_coefficients.ndim != 1 or len(self.electronic_coefficients) < 1:
            raise ValueError("variational SOC electronic coefficients must be nonempty.")
        if not all(
            np.all(np.isfinite(item))
            for item in (self.q, self.p, self.electronic_coefficients)
        ):
            raise ValueError("variational SOC state contains non-finite data.")
        if not np.isfinite(float(self.time_au)):
            raise ValueError("variational SOC time must be finite.")
        if abs(self.electronic_norm - 1.0) > float(tolerance):
            raise ValueError("variational SOC electronic coefficients must have unit norm.")
        return self

    def normalized(self):
        norm = self.electronic_norm
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError("cannot normalize a zero or non-finite electronic vector.")
        return CanonicalVariationalSOCStateV250(
            self.q,
            self.p,
            self.electronic_coefficients / np.sqrt(norm),
            self.time_au,
        ).validate()

    def as_dict(self):
        self.validate()
        return {
            "q": self.q.tolist(),
            "p": self.p.tolist(),
            "electronic_coefficients": _complex_pairs_v250(
                self.electronic_coefficients
            ),
            "time_au": float(self.time_au),
            "electronic_norm": self.electronic_norm,
        }


def _force_v250(point, coefficients, *, tolerance):
    coefficients = np.asarray(coefficients, dtype=complex)
    derivatives = np.asarray(point.hamiltonian_derivative_operator_q, dtype=complex)
    expectations = np.asarray(
        [np.vdot(coefficients, derivatives[a] @ coefficients) for a in range(point.nq)],
        dtype=complex,
    )
    if float(np.max(np.abs(expectations.imag))) > float(tolerance):
        raise ValueError("variational SOC force expectation has a non-negligible imaginary part.")
    return -expectations.real


def _energy_v250(state, mass_matrix, hamiltonian):
    kinetic = 0.5 * float(
        state.p @ np.linalg.solve(np.asarray(mass_matrix, dtype=float), state.p)
    )
    electronic = float(
        np.real(
            np.vdot(
                state.electronic_coefficients,
                np.asarray(hamiltonian, dtype=complex)
                @ state.electronic_coefficients,
            )
        )
    )
    return kinetic + electronic


@dataclass(frozen=True)
class SymmetricVariationalSOCStepV250:
    start: CanonicalVariationalSOCStateV250
    end: CanonicalVariationalSOCStateV250
    dt_au: float
    mass_matrix_au: np.ndarray
    mass_matrix_end_au: np.ndarray
    force_start: np.ndarray
    force_end: np.ndarray
    H_start: np.ndarray
    H_end: np.ndarray
    K_start: np.ndarray
    K_end: np.ndarray
    overlap_start_end: np.ndarray
    transport_end_to_start: np.ndarray
    singular_values: np.ndarray
    transport_policy: dict
    transport_metrics: dict
    energy_start_hartree: float
    energy_end_hartree: float

    def __post_init__(self):
        for name, dtype in (
            ("mass_matrix_au", float),
            ("mass_matrix_end_au", float),
            ("force_start", float),
            ("force_end", float),
            ("H_start", complex),
            ("H_end", complex),
            ("K_start", complex),
            ("K_end", complex),
            ("overlap_start_end", complex),
            ("transport_end_to_start", complex),
            ("singular_values", float),
        ):
            object.__setattr__(
                self, name, np.asarray(getattr(self, name), dtype=dtype).copy()
            )
        object.__setattr__(self, "transport_policy", dict(self.transport_policy))
        object.__setattr__(self, "transport_metrics", dict(self.transport_metrics))

    @property
    def energy_change_hartree(self):
        return float(self.energy_end_hartree - self.energy_start_hartree)

    @property
    def norm_change(self):
        return float(self.end.electronic_norm - self.start.electronic_norm)

    def validate(self, tolerance=2.0e-9):
        self.start.validate(tolerance=max(float(tolerance), 1.0e-10))
        self.end.validate(tolerance=max(float(tolerance), 1.0e-10))
        dt = float(self.dt_au)
        if not np.isfinite(dt) or dt == 0.0:
            raise ValueError("variational SOC step dt must be finite and nonzero.")
        nq = len(self.start.q)
        ns = len(self.start.electronic_coefficients)
        if self.end.q.shape != (nq,) or self.end.electronic_coefficients.shape != (ns,):
            raise ValueError("variational SOC state dimension changed within one step.")
        if self.mass_matrix_au.shape != (nq, nq):
            raise ValueError("variational SOC mass matrix has incompatible shape.")
        if self.mass_matrix_end_au.shape != (nq, nq):
            raise ValueError("variational SOC endpoint mass matrix has incompatible shape.")
        if not np.allclose(self.mass_matrix_au, self.mass_matrix_au.T, atol=tolerance):
            raise ValueError("variational SOC mass matrix is not symmetric.")
        if float(np.min(np.linalg.eigvalsh(self.mass_matrix_au))) <= 0.0:
            raise ValueError("variational SOC mass matrix is not positive definite.")
        if self.force_start.shape != (nq,) or self.force_end.shape != (nq,):
            raise ValueError("variational SOC forces have incompatible shape.")
        if self.H_start.shape != (ns, ns) or self.H_end.shape != (ns, ns):
            raise ValueError("variational SOC Hamiltonians have incompatible shape.")
        if self.K_start.shape != (nq, ns, ns) or self.K_end.shape != (nq, ns, ns):
            raise ValueError("variational SOC derivative operators have incompatible shape.")
        if self.overlap_start_end.shape != (ns, ns):
            raise ValueError("variational SOC overlap has incompatible shape.")
        if self.transport_end_to_start.shape != (ns, ns):
            raise ValueError("variational SOC polar transport has incompatible shape.")
        arrays = (
            self.mass_matrix_au,
            self.mass_matrix_end_au,
            self.force_start,
            self.force_end,
            self.H_start,
            self.H_end,
            self.K_start,
            self.K_end,
            self.overlap_start_end,
            self.transport_end_to_start,
            self.singular_values,
        )
        if any(not np.all(np.isfinite(item)) for item in arrays):
            raise ValueError("variational SOC step contains non-finite data.")
        if not np.allclose(
            self.mass_matrix_au,
            self.mass_matrix_end_au,
            rtol=1.0e-10,
            atol=1.0e-12,
        ):
            raise ValueError("velocity-Verlet step has a coordinate-dependent mass matrix.")
        for name, matrix in (("H_start", self.H_start), ("H_end", self.H_end)):
            if _scaled_frobenius_v250(matrix, matrix.conj().T) > tolerance:
                raise ValueError(f"{name} is not Hermitian.")
        for name, derivatives in (("K_start", self.K_start), ("K_end", self.K_end)):
            if any(
                _scaled_frobenius_v250(matrix, matrix.conj().T) > tolerance
                for matrix in derivatives
            ):
                raise ValueError(f"{name} contains a non-Hermitian component.")

        left, observed_singular, right_h = np.linalg.svd(
            self.overlap_start_end, full_matrices=False
        )
        expected_transport = left @ right_h
        if not np.allclose(
            self.singular_values, observed_singular, atol=tolerance, rtol=tolerance
        ):
            raise ValueError("stored singular values disagree with the raw overlap.")
        policy = self.transport_policy
        required_policy_keys = {
            "contraction_tolerance",
            "minimum_retained_singular_value",
            "maximum_condition_number",
            "maximum_principal_angle_radians",
            "transport_unitarity_tolerance",
            "polar_hermiticity_tolerance",
        }
        if set(policy) != required_policy_keys:
            raise ValueError("stored SVD-polar quality policy is incomplete.")
        minimum = float(np.min(observed_singular))
        maximum = float(np.max(observed_singular))
        condition = float(np.inf if minimum == 0.0 else maximum / minimum)
        principal_angle = float(np.arccos(np.clip(minimum, 0.0, 1.0)))
        if maximum - 1.0 > float(policy["contraction_tolerance"]):
            raise ValueError("raw variational SOC overlap has spectral expansion.")
        if minimum < float(policy["minimum_retained_singular_value"]):
            raise ValueError("raw variational SOC overlap loses the retained manifold.")
        if condition > float(policy["maximum_condition_number"]):
            raise ValueError("raw variational SOC overlap is ill-conditioned.")
        if principal_angle > float(policy["maximum_principal_angle_radians"]):
            raise ValueError("raw variational SOC overlap exceeds the principal-angle policy.")
        if _scaled_frobenius_v250(
            self.transport_end_to_start, expected_transport
        ) > tolerance:
            raise ValueError("stored polar transport disagrees with the overlap SVD.")
        identity = np.eye(ns, dtype=complex)
        if _scaled_frobenius_v250(
            self.transport_end_to_start.conj().T
            @ self.transport_end_to_start,
            identity,
        ) > tolerance:
            raise ValueError("stored polar transport is not unitary.")
        if self.transport_metrics.get("physically_consistent") is not True:
            raise ValueError("variational SOC raw overlap is not physically consistent.")
        if self.transport_metrics.get("trajectory_ready") is not True:
            raise ValueError("variational SOC raw overlap is not trajectory quality.")
        metric_singular_values = np.asarray(
            self.transport_metrics.get("singular_values", ()), dtype=float
        )
        if not np.allclose(
            metric_singular_values,
            observed_singular,
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("transport metrics are not bound to the overlap singular values.")
        for key, observed in (
            ("minimum_singular_value", minimum),
            ("maximum_singular_value", maximum),
            ("condition_number", condition),
            ("maximum_principal_angle_radians", principal_angle),
        ):
            if abs(float(self.transport_metrics.get(key, np.inf)) - observed) > tolerance:
                raise ValueError(f"transport metric {key} disagrees with the raw overlap.")

        expected_force_start = -np.asarray(
            [
                np.real(
                    np.vdot(
                        self.start.electronic_coefficients,
                        self.K_start[a] @ self.start.electronic_coefficients,
                    )
                )
                for a in range(nq)
            ]
        )
        expected_force_end = -np.asarray(
            [
                np.real(
                    np.vdot(
                        self.end.electronic_coefficients,
                        self.K_end[a] @ self.end.electronic_coefficients,
                    )
                )
                for a in range(nq)
            ]
        )
        if not np.allclose(self.force_start, expected_force_start, atol=tolerance):
            raise ValueError("stored start force disagrees with the variational state.")
        if not np.allclose(self.force_end, expected_force_end, atol=tolerance):
            raise ValueError("stored endpoint force disagrees with the variational state.")

        p_half = self.start.p + 0.5 * dt * self.force_start
        expected_q = self.start.q + dt * np.linalg.solve(self.mass_matrix_au, p_half)
        expected_p = p_half + 0.5 * dt * self.force_end
        if not np.allclose(self.end.q, expected_q, atol=tolerance, rtol=tolerance):
            raise ValueError("stored endpoint coordinate disagrees with velocity-Verlet.")
        if not np.allclose(self.end.p, expected_p, atol=tolerance, rtol=tolerance):
            raise ValueError("stored endpoint momentum disagrees with velocity-Verlet.")
        expected_c = (
            hermitian_exponential(self.H_end, 0.5 * dt)
            @ self.transport_end_to_start.conj().T
            @ hermitian_exponential(self.H_start, 0.5 * dt)
            @ self.start.electronic_coefficients
        )
        if not np.allclose(
            self.end.electronic_coefficients,
            expected_c,
            atol=tolerance,
            rtol=tolerance,
        ):
            raise ValueError("stored electronic endpoint disagrees with Strang-polar propagation.")
        if abs(self.end.time_au - (self.start.time_au + dt)) > tolerance:
            raise ValueError("stored endpoint time disagrees with dt.")

        expected_energy_start = _energy_v250(
            self.start, self.mass_matrix_au, self.H_start
        )
        expected_energy_end = _energy_v250(self.end, self.mass_matrix_au, self.H_end)
        if abs(float(self.energy_start_hartree) - expected_energy_start) > tolerance:
            raise ValueError("stored start energy is inconsistent.")
        if abs(float(self.energy_end_hartree) - expected_energy_end) > tolerance:
            raise ValueError("stored endpoint energy is inconsistent.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "start": self.start.as_dict(),
            "end": self.end.as_dict(),
            "dt_au": float(self.dt_au),
            "mass_matrix_au": self.mass_matrix_au.tolist(),
            "mass_matrix_end_au": self.mass_matrix_end_au.tolist(),
            "force_start": self.force_start.tolist(),
            "force_end": self.force_end.tolist(),
            "H_start": _complex_pairs_v250(self.H_start),
            "H_end": _complex_pairs_v250(self.H_end),
            "K_start": _complex_pairs_v250(self.K_start),
            "K_end": _complex_pairs_v250(self.K_end),
            "overlap_start_end": _complex_pairs_v250(self.overlap_start_end),
            "transport_end_to_start": _complex_pairs_v250(
                self.transport_end_to_start
            ),
            "singular_values": self.singular_values.tolist(),
            "transport_policy": self.transport_policy,
            "transport_metrics": self.transport_metrics,
            "energy_start_hartree": float(self.energy_start_hartree),
            "energy_end_hartree": float(self.energy_end_hartree),
            "energy_change_hartree": self.energy_change_hartree,
            "norm_change": self.norm_change,
        }


def symmetric_variational_soc_step_v250(
    state,
    provider,
    dt_au,
    *,
    settings=VariationalSOCIntegratorSettingsV250(),
):
    """Advance one self-adjoint restricted-TDVP step.

    A signed ``dt_au`` is accepted so forward/backward reversibility can be tested
    without defining a distinct reverse algorithm.
    """

    settings = settings.validate()
    state = state.validate(tolerance=settings.structural_tolerance)
    dt = float(dt_au)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError("variational SOC step dt must be finite and nonzero.")
    for method_name in ("evaluate_snapshot", "snapshot_overlap"):
        if not callable(getattr(provider, method_name, None)):
            raise TypeError(
                "variational SOC provider must implement evaluate_snapshot and snapshot_overlap."
            )

    snapshot_start = provider.evaluate_snapshot(state.q)
    if not hasattr(snapshot_start, "point"):
        raise TypeError(
            "v0.25.0 requires an ElectronicOperatorSnapshot with full H, K, D, and mass."
        )
    snapshot_start = snapshot_start.validate(
        atol=settings.structural_tolerance,
        isometry_atol=settings.structural_tolerance,
    )
    point_start = snapshot_start.point
    if not np.allclose(point_start.q, state.q, atol=settings.structural_tolerance):
        raise ValueError("provider start snapshot geometry disagrees with the state.")
    if point_start.nstate != len(state.electronic_coefficients):
        raise ValueError("provider electronic dimension disagrees with the state.")

    mass_start = np.asarray(point_start.mass_matrix_q_au, dtype=float)
    force_start = _force_v250(
        point_start,
        state.electronic_coefficients,
        tolerance=settings.structural_tolerance,
    )
    p_half = state.p + 0.5 * dt * force_start
    q_end = state.q + dt * np.linalg.solve(mass_start, p_half)

    snapshot_end = provider.evaluate_snapshot(q_end)
    if not hasattr(snapshot_end, "point"):
        raise TypeError("provider endpoint is not an ElectronicOperatorSnapshot.")
    snapshot_end = snapshot_end.validate(
        atol=settings.structural_tolerance,
        isometry_atol=settings.structural_tolerance,
    )
    point_end = snapshot_end.point
    if not np.allclose(point_end.q, q_end, atol=settings.structural_tolerance):
        raise ValueError("provider endpoint snapshot geometry disagrees with Verlet drift.")
    mass_end = np.asarray(point_end.mass_matrix_q_au, dtype=float)
    if not np.allclose(
        mass_start,
        mass_end,
        rtol=settings.mass_relative_tolerance,
        atol=settings.mass_absolute_tolerance,
    ):
        raise ValueError(
            "v0.25.0 velocity-Verlet requires a constant generalized mass matrix."
        )

    overlap = np.asarray(
        provider.snapshot_overlap(snapshot_start, snapshot_end), dtype=complex
    )
    transport = certified_transport_from_overlap_v233(
        overlap, policy=settings.overlap_policy
    )
    W = transport.right_to_left_transport
    H_start = np.asarray(point_start.H, dtype=complex)
    H_end = np.asarray(point_end.H, dtype=complex)
    c_end = (
        hermitian_exponential(H_end, 0.5 * dt)
        @ W.conj().T
        @ hermitian_exponential(H_start, 0.5 * dt)
        @ state.electronic_coefficients
    )
    provisional_end = CanonicalVariationalSOCStateV250(
        q_end,
        p_half,
        c_end,
        state.time_au + dt,
    ).validate(tolerance=settings.step_binding_tolerance)
    force_end = _force_v250(
        point_end,
        provisional_end.electronic_coefficients,
        tolerance=settings.structural_tolerance,
    )
    p_end = p_half + 0.5 * dt * force_end
    end = CanonicalVariationalSOCStateV250(
        q_end,
        p_end,
        c_end,
        state.time_au + dt,
    ).validate(tolerance=settings.step_binding_tolerance)
    return SymmetricVariationalSOCStepV250(
        start=state,
        end=end,
        dt_au=dt,
        mass_matrix_au=mass_start,
        mass_matrix_end_au=mass_end,
        force_start=force_start,
        force_end=force_end,
        H_start=H_start,
        H_end=H_end,
        K_start=point_start.hamiltonian_derivative_operator_q,
        K_end=point_end.hamiltonian_derivative_operator_q,
        overlap_start_end=overlap,
        transport_end_to_start=W,
        singular_values=transport.singular_values,
        transport_policy=settings.overlap_policy.as_dict(),
        transport_metrics=transport.as_dict(),
        energy_start_hartree=_energy_v250(state, mass_start, H_start),
        energy_end_hartree=_energy_v250(end, mass_start, H_end),
    ).validate(tolerance=settings.step_binding_tolerance)


V250_TRAJECTORY_CLAIMS = {
    "restricted_single_packet_tdvp_validated": True,
    "symmetric_strang_verlet_coupling_validated": True,
    "svd_computed_polar_transport_validated": True,
    "complete_spinor_soc_propagation_validated": True,
    "coordinate_dependent_complex_gauge_covariance_validated": True,
    "full_multi_gaussian_tdvp_validated": False,
    "adaptive_gaussian_width_tdvp_validated": False,
    "plain_verlet_for_general_tdvp_validated": False,
    "coordinate_dependent_mass_verlet_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
    "general_ab_initio_soc_dynamics_accuracy_validated": False,
}


@dataclass(frozen=True)
class SymmetricVariationalSOCTrajectoryV250:
    initial_state: CanonicalVariationalSOCStateV250
    final_state: CanonicalVariationalSOCStateV250
    steps: tuple
    settings: VariationalSOCIntegratorSettingsV250
    claims: dict

    def __post_init__(self):
        object.__setattr__(self, "steps", tuple(self.steps))
        object.__setattr__(self, "claims", dict(self.claims))

    @property
    def maximum_norm_drift(self):
        states = [self.initial_state] + [step.end for step in self.steps]
        return float(max(abs(state.electronic_norm - 1.0) for state in states))

    @property
    def maximum_absolute_energy_drift_hartree(self):
        if not self.steps:
            return 0.0
        reference = self.steps[0].energy_start_hartree
        return float(
            max(abs(step.energy_end_hartree - reference) for step in self.steps)
        )

    def validate(self):
        settings = self.settings.validate()
        self.initial_state.validate(tolerance=settings.step_binding_tolerance)
        self.final_state.validate(tolerance=settings.step_binding_tolerance)
        if any(type(value) is not bool for value in self.claims.values()):
            raise TypeError("every v0.25.0 trajectory claim must be a native Boolean.")
        if self.claims != V250_TRAJECTORY_CLAIMS:
            raise ValueError("v0.25.0 trajectory claims differ from the frozen boundary.")
        previous = self.initial_state
        for step in self.steps:
            step.validate(tolerance=settings.step_binding_tolerance)
            if _state_distance_v250(step.start, previous) > settings.step_binding_tolerance:
                raise ValueError("v0.25.0 trajectory step chain is discontinuous.")
            previous = step.end
        if _state_distance_v250(previous, self.final_state) > settings.step_binding_tolerance:
            raise ValueError("v0.25.0 final state is not bound to the last step.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": VARIATIONAL_SOC_SCHEMA_V250,
            "ansatz": RESTRICTED_TDVP_ANSATZ_V250,
            "general_tdvp_integrator_recommendation": GENERAL_TDVP_INTEGRATOR_V250,
            "settings": self.settings.as_dict(),
            "initial_state": self.initial_state.as_dict(),
            "final_state": self.final_state.as_dict(),
            "steps": [step.as_dict() for step in self.steps],
            "maximum_norm_drift": self.maximum_norm_drift,
            "maximum_absolute_energy_drift_hartree": (
                self.maximum_absolute_energy_drift_hartree
            ),
            "claims": dict(self.claims),
        }

    def fingerprint(self):
        return _sha256_v250(self.as_dict())


def run_symmetric_variational_soc_dynamics_v250(
    initial_state,
    provider,
    *,
    dt_au,
    steps,
    settings=VariationalSOCIntegratorSettingsV250(),
):
    """Run the validated forward v0.25.0 restricted-TDVP trajectory."""

    settings = settings.validate()
    dt = float(dt_au)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("forward variational SOC dt must be finite and positive.")
    if int(steps) != steps or int(steps) < 0:
        raise ValueError("variational SOC steps must be a nonnegative integer.")
    state = initial_state.normalized()
    initial = state
    receipts = []
    for _ in range(int(steps)):
        receipt = symmetric_variational_soc_step_v250(
            state, provider, dt, settings=settings
        )
        receipts.append(receipt)
        state = receipt.end
    return SymmetricVariationalSOCTrajectoryV250(
        initial_state=initial,
        final_state=state,
        steps=tuple(receipts),
        settings=settings,
        claims=V250_TRAJECTORY_CLAIMS,
    ).validate()


def reverse_variational_soc_trajectory_v250(trajectory, provider):
    """Apply the exact signed-step adjoint sequence for reversibility audits."""

    trajectory = trajectory.validate()
    state = trajectory.final_state
    reverse_steps = []
    for forward_step in reversed(trajectory.steps):
        reverse_step = symmetric_variational_soc_step_v250(
            state,
            provider,
            -forward_step.dt_au,
            settings=trajectory.settings,
        )
        reverse_steps.append(reverse_step)
        state = reverse_step.end
    return SymmetricVariationalSOCTrajectoryV250(
        initial_state=trajectory.final_state,
        final_state=state,
        steps=tuple(reverse_steps),
        settings=trajectory.settings,
        claims=V250_TRAJECTORY_CLAIMS,
    ).validate()
