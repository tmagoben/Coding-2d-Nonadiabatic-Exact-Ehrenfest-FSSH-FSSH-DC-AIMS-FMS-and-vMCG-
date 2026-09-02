"""Controlled adaptive-basis lifecycle for v0.25.3.

This layer leaves the v0.25.2 adaptive-width TDVP vector field unchanged.  At a
configured checkpoint it permits at most one auditable topology event:

* merge one nearly redundant packet into a retained survivor,
* prune one low-population packet when an SVD projection proves it dispensable, or
* spawn one analytically scored residual direction.

Every accepted event is a full-SVD variational projection at fixed time.  Projection
loss, fidelity, norm, energy jump, basis rank, condition number, stable packet IDs,
and packet ages are bound into a self-validating receipt.  The released scope remains
one nuclear coordinate, scalar positive widths with chirps, a fixed complete
electronic frame, and a Hermitian quadratic matrix Hamiltonian.
"""

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.optimize import root

from .adaptive_multigaussian_tdvp_v252 import (
    AdaptiveVariationalSettingsV252,
    QuadraticSpinHamiltonianV252,
    ThawedGaussianSpinorStateV252,
    _complex_pairs_v252,
    _cross_moments_v252,
    _integrate_polynomial_v252,
    _kinetic_polynomial_v252,
    _scaled_norm_v252,
    _sha256_v252,
    _tangent_terms_v252,
    _validate_width_domain_v252,
    adaptive_implicit_midpoint_tdvp_step_v252,
    adaptive_variational_energy_v252,
    build_adaptive_variational_metric_system_v252,
    pack_adaptive_variational_parameters_v252,
    quadratic_spin_hamiltonian_from_provider_v252,
    state_from_adaptive_variational_parameters_v252,
)
from .multigaussian_tdvp_v251 import MetricSolveReceiptV251, solve_variational_metric_v251


CONTROLLED_BASIS_SCHEMA_V253 = "gnd-controlled-basis-trajectory-v0.25.3"
CONTROLLED_BASIS_EVENT_SCHEMA_V253 = "gnd-controlled-basis-event-v0.25.3"
SPAWN_SCORE_V253 = (
    "norm of <orthogonalized candidate electronic spinors|dPsi/dt+iHPsi>"
)
PROJECTION_POLICY_V253 = (
    "full-SVD least-squares projection with explicit loss, rank, conditioning, "
    "fidelity, norm, and energy-jump gates"
)
EVENT_ORDER_V253 = (
    "at most one event per checkpoint: merge, then prune, then residual spawn"
)
POTENTIAL_CONTRACT_V253 = (
    "fixed-frame one-dimensional Hermitian quadratic complete-spinor Hamiltonian"
)


def _positive_integer_v253(value, name, *, allow_zero=False):
    lower = 0 if allow_zero else 1
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        or int(value) < lower
    ):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be a {qualifier} integer.")
    return int(value)


def _basis_arrays_v253(q, p, widths, chirps):
    arrays = tuple(np.asarray(value, dtype=float).copy() for value in (q, p, widths, chirps))
    if arrays[0].ndim != 1 or len(arrays[0]) < 1:
        raise ValueError("target Gaussian basis must be a nonempty vector.")
    if any(value.shape != arrays[0].shape for value in arrays[1:]):
        raise ValueError("target Gaussian basis arrays have different shapes.")
    if not all(np.all(np.isfinite(value)) for value in arrays):
        raise ValueError("target Gaussian basis contains non-finite values.")
    if np.any(arrays[2] <= 0.0):
        raise ValueError("target Gaussian widths must be positive.")
    return arrays


def _state_distance_v253(left, right):
    if left.ngaussian != right.ngaussian or left.nstate != right.nstate:
        return float("inf")
    return max(
        _scaled_norm_v252(
            pack_adaptive_variational_parameters_v252(left),
            pack_adaptive_variational_parameters_v252(right),
        ),
        abs(float(left.time_au) - float(right.time_au)),
    )


def _validate_packet_metadata_v253(packet_ids, packet_ages, count):
    packet_ids = tuple(packet_ids)
    packet_ages = tuple(packet_ages)
    if len(packet_ids) != count or len(packet_ages) != count:
        raise ValueError("packet metadata length disagrees with Gaussian count.")
    if any(not isinstance(value, str) or not value for value in packet_ids):
        raise ValueError("packet IDs must be nonempty strings.")
    if len(set(packet_ids)) != len(packet_ids):
        raise ValueError("packet IDs must be unique.")
    ages = tuple(_positive_integer_v253(value, "packet age", allow_zero=True) for value in packet_ages)
    return packet_ids, ages


@dataclass(frozen=True)
class ControlledBasisSettingsV253:
    """Frozen topology algorithms and conservative admission gates."""

    tdvp_settings: AdaptiveVariationalSettingsV252 = field(
        default_factory=AdaptiveVariationalSettingsV252
    )
    adaptation_interval_steps: int = 1
    minimum_packet_count: int = 1
    maximum_packet_count: int = 8
    minimum_packet_age_steps: int = 2
    maximum_activation_age_steps: int = 64
    position_displacement_sigma: float = 2.0
    momentum_displacement_sigma: float = 2.0
    spawn_residual_capture_threshold: float = 1.0e-5
    minimum_spawn_novelty: float = 1.0e-4
    projection_relative_cutoff: float = 1.0e-11
    projection_absolute_cutoff: float = 1.0e-13
    maximum_basis_condition_number: float = 1.0e8
    projection_linear_residual_tolerance: float = 2.0e-10
    maximum_spawn_projection_loss: float = 2.0e-10
    maximum_prune_coefficient_population: float = 1.0e-8
    maximum_prune_projection_loss: float = 1.0e-8
    minimum_merge_overlap: float = 0.999
    maximum_merge_projection_loss: float = 1.0e-8
    maximum_event_energy_jump_hartree: float = 1.0e-6
    shape_activation_population: float = 1.0e-6
    structural_tolerance: float = 2.0e-9
    spawning: bool = True
    pruning: bool = True
    merging: bool = True
    one_event_per_checkpoint: bool = True
    residual_score: str = SPAWN_SCORE_V253
    projection_policy: str = PROJECTION_POLICY_V253
    event_order: str = EVENT_ORDER_V253
    multidimensional_nuclear_motion: bool = False
    full_width_matrices: bool = False
    coordinate_dependent_electronic_frame: bool = False
    real_molecular_soc_provider: bool = False
    general_aims_branching: bool = False

    def validate(self):
        self.tdvp_settings.validate()
        interval = _positive_integer_v253(
            self.adaptation_interval_steps, "adaptation_interval_steps"
        )
        minimum = _positive_integer_v253(self.minimum_packet_count, "minimum_packet_count")
        maximum = _positive_integer_v253(self.maximum_packet_count, "maximum_packet_count")
        _positive_integer_v253(
            self.minimum_packet_age_steps, "minimum_packet_age_steps", allow_zero=True
        )
        maximum_activation_age = _positive_integer_v253(
            self.maximum_activation_age_steps, "maximum_activation_age_steps"
        )
        if maximum_activation_age <= int(self.minimum_packet_age_steps):
            raise ValueError(
                "maximum_activation_age_steps must exceed minimum_packet_age_steps."
            )
        if minimum > maximum:
            raise ValueError("minimum_packet_count cannot exceed maximum_packet_count.")
        if interval < 1:
            raise ValueError("adaptation interval must be positive.")
        for name in (
            "position_displacement_sigma",
            "momentum_displacement_sigma",
            "spawn_residual_capture_threshold",
            "minimum_spawn_novelty",
            "projection_relative_cutoff",
            "projection_absolute_cutoff",
            "maximum_basis_condition_number",
            "projection_linear_residual_tolerance",
            "maximum_spawn_projection_loss",
            "maximum_prune_coefficient_population",
            "maximum_prune_projection_loss",
            "minimum_merge_overlap",
            "maximum_merge_projection_loss",
            "maximum_event_energy_jump_hartree",
            "shape_activation_population",
            "structural_tolerance",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if float(self.projection_relative_cutoff) >= 1.0:
            raise ValueError("projection_relative_cutoff must be smaller than one.")
        if float(self.minimum_spawn_novelty) >= 1.0:
            raise ValueError("minimum_spawn_novelty must be smaller than one.")
        if float(self.minimum_merge_overlap) >= 1.0:
            raise ValueError("minimum_merge_overlap must be smaller than one.")
        if float(self.maximum_basis_condition_number) < 1.0:
            raise ValueError("maximum_basis_condition_number must be at least one.")
        for name in (
            "spawning",
            "pruning",
            "merging",
            "one_event_per_checkpoint",
            "multidimensional_nuclear_motion",
            "full_width_matrices",
            "coordinate_dependent_electronic_frame",
            "real_molecular_soc_provider",
            "general_aims_branching",
        ):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
        if not self.one_event_per_checkpoint:
            raise ValueError("v0.25.3 freezes one topology event per checkpoint.")
        if not (self.spawning and self.pruning and self.merging):
            raise ValueError("v0.25.3 freezes spawn, prune, and merge procedures as enabled.")
        closed = {
            "multidimensional_nuclear_motion": self.multidimensional_nuclear_motion,
            "full_width_matrices": self.full_width_matrices,
            "coordinate_dependent_electronic_frame": self.coordinate_dependent_electronic_frame,
            "real_molecular_soc_provider": self.real_molecular_soc_provider,
            "general_aims_branching": self.general_aims_branching,
        }
        requested = sorted(name for name, enabled in closed.items() if enabled)
        if requested:
            raise ValueError("v0.25.3 does not admit: " + ", ".join(requested) + ".")
        if self.residual_score != SPAWN_SCORE_V253:
            raise ValueError("the v0.25.3 residual score is frozen.")
        if self.projection_policy != PROJECTION_POLICY_V253:
            raise ValueError("the v0.25.3 projection policy is frozen.")
        if self.event_order != EVENT_ORDER_V253:
            raise ValueError("the v0.25.3 event order is frozen.")
        return self

    def as_dict(self):
        self.validate()
        result = asdict(self)
        result["tdvp_settings"] = self.tdvp_settings.as_dict()
        return result


@dataclass(frozen=True)
class SpawnCandidateV253:
    q: float
    p: float
    width: float
    chirp: float
    parent_index: int
    displacement_kind: str

    def validate(self, *, packet_count=None):
        for name in ("q", "p", "width", "chirp"):
            if not np.isfinite(float(getattr(self, name))):
                raise ValueError(f"spawn candidate {name} must be finite.")
        if float(self.width) <= 0.0:
            raise ValueError("spawn candidate width must be positive.")
        parent = _positive_integer_v253(self.parent_index, "parent_index", allow_zero=True)
        if packet_count is not None and parent >= int(packet_count):
            raise ValueError("spawn candidate parent index is outside the basis.")
        if self.displacement_kind not in (
            "position-minus",
            "position-plus",
            "momentum-minus",
            "momentum-plus",
            "external",
        ):
            raise ValueError("spawn candidate displacement kind is unsupported.")
        return self

    def canonical_key(self):
        self.validate()
        return (
            float(self.q),
            float(self.p),
            float(self.width),
            float(self.chirp),
            self.displacement_kind,
        )

    def as_dict(self):
        self.validate()
        return {
            "q": float(self.q),
            "p": float(self.p),
            "width": float(self.width),
            "chirp": float(self.chirp),
            "parent_index": int(self.parent_index),
            "displacement_kind": self.displacement_kind,
        }


def generate_spawn_candidates_v253(state, *, settings=ControlledBasisSettingsV253()):
    state = state.validate()
    settings = settings.validate()
    candidates = []
    for index in range(state.ngaussian):
        sigma_q = 1.0 / np.sqrt(2.0 * state.widths[index])
        sigma_p = np.sqrt(state.widths[index] / 2.0)
        base = dict(
            width=float(state.widths[index]),
            chirp=float(state.chirps[index]),
            parent_index=index,
        )
        candidates.extend(
            (
                SpawnCandidateV253(
                    q=float(state.q[index] - settings.position_displacement_sigma * sigma_q),
                    p=float(state.p[index]),
                    displacement_kind="position-minus",
                    **base,
                ),
                SpawnCandidateV253(
                    q=float(state.q[index] + settings.position_displacement_sigma * sigma_q),
                    p=float(state.p[index]),
                    displacement_kind="position-plus",
                    **base,
                ),
                SpawnCandidateV253(
                    q=float(state.q[index]),
                    p=float(state.p[index] - settings.momentum_displacement_sigma * sigma_p),
                    displacement_kind="momentum-minus",
                    **base,
                ),
                SpawnCandidateV253(
                    q=float(state.q[index]),
                    p=float(state.p[index] + settings.momentum_displacement_sigma * sigma_p),
                    displacement_kind="momentum-plus",
                    **base,
                ),
            )
        )
    return tuple(candidate.validate(packet_count=state.ngaussian) for candidate in candidates)


def _cross_overlap_v253(bra_q, bra_p, bra_widths, bra_chirps, ket_state):
    count = len(bra_q)
    result = np.zeros((count, ket_state.ngaussian), dtype=complex)
    for i in range(count):
        for j in range(ket_state.ngaussian):
            result[i, j] = _cross_moments_v252(
                bra_q[i], bra_p[i], bra_widths[i], bra_chirps[i],
                ket_state.q[j], ket_state.p[j], ket_state.widths[j], ket_state.chirps[j],
                maximum_order=0,
            )[0]
    return result


def _nuclear_overlap_from_arrays_v253(q, p, widths, chirps):
    count = len(q)
    result = np.zeros((count, count), dtype=complex)
    for i in range(count):
        for j in range(count):
            result[i, j] = _cross_moments_v252(
                q[i], p[i], widths[i], chirps[i],
                q[j], p[j], widths[j], chirps[j],
                maximum_order=0,
            )[0]
    return result


def _svd_basis_solve_v253(matrix, rhs, settings):
    matrix = np.asarray(matrix, dtype=complex)
    rhs = np.asarray(rhs, dtype=complex)
    u, singular_values, vh = np.linalg.svd(matrix, full_matrices=True)
    largest = float(singular_values[0]) if len(singular_values) else 0.0
    cutoff = max(
        float(settings.projection_absolute_cutoff),
        float(settings.projection_relative_cutoff) * largest,
    )
    retained = singular_values > cutoff
    rank = int(np.count_nonzero(retained))
    if rank != matrix.shape[0]:
        raise ValueError("target Gaussian basis is rank deficient under the projection cutoff.")
    condition = float(singular_values[0] / singular_values[-1])
    if condition > settings.maximum_basis_condition_number:
        raise ValueError("target Gaussian basis exceeds the conditioning gate.")
    solution = vh.conj().T @ ((u.conj().T @ rhs) / singular_values[:, None])
    residual = matrix @ solution - rhs
    residual_relative = float(
        np.linalg.norm(residual) / max(float(np.linalg.norm(rhs)), 1.0e-30)
    )
    if residual_relative > settings.projection_linear_residual_tolerance:
        raise ValueError("basis projection linear residual exceeds tolerance.")
    return solution, singular_values, cutoff, condition, residual, residual_relative


def _projection_data_v253(source, q, p, widths, chirps, model, settings):
    source = source.validate()
    model = model.validate()
    settings = settings.validate()
    q, p, widths, chirps = _basis_arrays_v253(q, p, widths, chirps)
    if source.nstate != model.nstate:
        raise ValueError("projection source/model electronic dimensions differ.")
    if np.min(widths) < settings.tdvp_settings.minimum_width or np.max(widths) > settings.tdvp_settings.maximum_width:
        raise ValueError("projection target width is outside the TDVP domain.")
    if np.max(np.abs(chirps)) > settings.tdvp_settings.maximum_absolute_chirp:
        raise ValueError("projection target chirp is outside the TDVP domain.")
    overlap = _nuclear_overlap_from_arrays_v253(q, p, widths, chirps)
    cross = _cross_overlap_v253(q, p, widths, chirps, source)
    rhs = cross @ source.coefficients
    raw_coefficients, singular_values, cutoff, condition, linear_residual, linear_relative = (
        _svd_basis_solve_v253(overlap, rhs, settings)
    )
    raw_state = ThawedGaussianSpinorStateV252(
        q=q, p=p, widths=widths, chirps=chirps,
        coefficients=raw_coefficients, time_au=source.time_au,
    ).validate(require_normalized=False)
    source_norm = source.generalized_norm
    raw_norm = raw_state.generalized_norm
    cross_inner = np.vdot(raw_coefficients, rhs)
    residual_squared = float(
        source_norm - 2.0 * np.real(cross_inner) + raw_norm
    )
    tolerance = settings.structural_tolerance * max(source_norm, 1.0)
    if residual_squared < -tolerance:
        raise ValueError("basis projection produced a negative residual norm.")
    residual_squared = max(residual_squared, 0.0)
    relative_loss = residual_squared / source_norm
    projected = raw_state.normalized()
    normalized_inner = cross_inner / np.sqrt(source_norm * raw_norm)
    fidelity = float(abs(normalized_inner) ** 2)
    source_energy = adaptive_variational_energy_v252(source, model)
    target_energy = adaptive_variational_energy_v252(projected, model)
    return {
        "overlap": overlap,
        "cross": cross,
        "rhs": rhs,
        "raw_coefficients": raw_coefficients,
        "singular_values": singular_values,
        "cutoff": cutoff,
        "condition": condition,
        "linear_residual": linear_residual,
        "linear_relative": linear_relative,
        "source_norm": source_norm,
        "raw_norm": raw_norm,
        "relative_loss": relative_loss,
        "fidelity": fidelity,
        "projected": projected,
        "source_energy": source_energy,
        "target_energy": target_energy,
    }


@dataclass(frozen=True)
class BasisProjectionReceiptV253:
    event_kind: str
    source: ThawedGaussianSpinorStateV252
    projected: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV252
    settings: ControlledBasisSettingsV253
    target_overlap: np.ndarray
    source_target_cross_overlap: np.ndarray
    projection_rhs: np.ndarray
    raw_projected_coefficients: np.ndarray
    singular_values: np.ndarray
    cutoff: float
    retained_condition_number: float
    linear_residual: np.ndarray
    linear_residual_relative: float
    source_norm: float
    raw_projected_norm: float
    relative_projection_loss: float
    normalized_fidelity: float
    source_energy_hartree: float
    projected_energy_hartree: float

    def __post_init__(self):
        for name in (
            "target_overlap", "source_target_cross_overlap", "projection_rhs",
            "raw_projected_coefficients", "linear_residual",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=complex).copy())
        object.__setattr__(self, "singular_values", np.asarray(self.singular_values, dtype=float).copy())

    @property
    def energy_jump_hartree(self):
        return float(self.projected_energy_hartree - self.source_energy_hartree)

    def validate(self):
        if self.event_kind not in ("spawn", "prune", "merge"):
            raise ValueError("basis projection event kind is invalid.")
        expected = _projection_data_v253(
            self.source, self.projected.q, self.projected.p,
            self.projected.widths, self.projected.chirps, self.model, self.settings,
        )
        tolerance = self.settings.structural_tolerance
        array_fields = {
            "target_overlap": "overlap",
            "source_target_cross_overlap": "cross",
            "projection_rhs": "rhs",
            "raw_projected_coefficients": "raw_coefficients",
            "singular_values": "singular_values",
            "linear_residual": "linear_residual",
        }
        for stored_name, expected_name in array_fields.items():
            if _scaled_norm_v252(getattr(self, stored_name), expected[expected_name]) > tolerance:
                raise ValueError(f"stored basis projection {stored_name} is inconsistent.")
        if _state_distance_v253(self.projected, expected["projected"]) > tolerance:
            raise ValueError("stored projected state is inconsistent.")
        scalar_fields = {
            "cutoff": "cutoff",
            "retained_condition_number": "condition",
            "linear_residual_relative": "linear_relative",
            "source_norm": "source_norm",
            "raw_projected_norm": "raw_norm",
            "relative_projection_loss": "relative_loss",
            "normalized_fidelity": "fidelity",
            "source_energy_hartree": "source_energy",
            "projected_energy_hartree": "target_energy",
        }
        for stored_name, expected_name in scalar_fields.items():
            if abs(float(getattr(self, stored_name)) - float(expected[expected_name])) > tolerance:
                raise ValueError(f"stored basis projection {stored_name} is inconsistent.")
        loss_limit = {
            "spawn": self.settings.maximum_spawn_projection_loss,
            "prune": self.settings.maximum_prune_projection_loss,
            "merge": self.settings.maximum_merge_projection_loss,
        }[self.event_kind]
        if self.relative_projection_loss > loss_limit:
            raise ValueError(f"{self.event_kind} projection loss exceeds its release gate.")
        if abs(self.energy_jump_hartree) > self.settings.maximum_event_energy_jump_hartree:
            raise ValueError("basis projection energy jump exceeds the release gate.")
        if abs(self.projected.generalized_norm - 1.0) > tolerance:
            raise ValueError("projected adaptive state is not normalized.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "event_kind": self.event_kind,
            "source": self.source.as_dict(),
            "projected": self.projected.as_dict(),
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "target_overlap": _complex_pairs_v252(self.target_overlap),
            "source_target_cross_overlap": _complex_pairs_v252(self.source_target_cross_overlap),
            "projection_rhs": _complex_pairs_v252(self.projection_rhs),
            "raw_projected_coefficients": _complex_pairs_v252(self.raw_projected_coefficients),
            "singular_values": self.singular_values.tolist(),
            "cutoff": float(self.cutoff),
            "retained_condition_number": float(self.retained_condition_number),
            "linear_residual": _complex_pairs_v252(self.linear_residual),
            "linear_residual_relative": float(self.linear_residual_relative),
            "source_norm": float(self.source_norm),
            "raw_projected_norm": float(self.raw_projected_norm),
            "relative_projection_loss": float(self.relative_projection_loss),
            "normalized_fidelity": float(self.normalized_fidelity),
            "source_energy_hartree": float(self.source_energy_hartree),
            "projected_energy_hartree": float(self.projected_energy_hartree),
            "energy_jump_hartree": self.energy_jump_hartree,
        }


def project_adaptive_state_v253(
    source, model, *, q, p, widths, chirps, event_kind,
    settings=ControlledBasisSettingsV253(),
):
    data = _projection_data_v253(source, q, p, widths, chirps, model, settings)
    return BasisProjectionReceiptV253(
        event_kind=event_kind,
        source=source,
        projected=data["projected"],
        model=model,
        settings=settings,
        target_overlap=data["overlap"],
        source_target_cross_overlap=data["cross"],
        projection_rhs=data["rhs"],
        raw_projected_coefficients=data["raw_coefficients"],
        singular_values=data["singular_values"],
        cutoff=data["cutoff"],
        retained_condition_number=data["condition"],
        linear_residual=data["linear_residual"],
        linear_residual_relative=data["linear_relative"],
        source_norm=data["source_norm"],
        raw_projected_norm=data["raw_norm"],
        relative_projection_loss=data["relative_loss"],
        normalized_fidelity=data["fidelity"],
        source_energy_hartree=data["source_energy"],
        projected_energy_hartree=data["target_energy"],
    ).validate()


def _residual_coupling_for_geometry_v253(state, model, velocity, geometry):
    """Return <g_c,a|dot(Psi)+i H Psi> for every electronic component."""

    q_c, p_c, width_c, chirp_c = geometry
    tangent_overlaps = np.zeros((state.nstate, state.parameter_count), dtype=complex)
    for mu, (packet, electronic_vector, polynomial) in enumerate(_tangent_terms_v252(state)):
        moments = _cross_moments_v252(
            q_c, p_c, width_c, chirp_c,
            state.q[packet], state.p[packet], state.widths[packet], state.chirps[packet],
            maximum_order=len(polynomial) - 1,
        )
        tangent_overlaps[:, mu] = electronic_vector * _integrate_polynomial_v252(
            polynomial, moments
        )
    dot_overlap = tangent_overlaps @ velocity
    h_overlap = np.zeros(state.nstate, dtype=complex)
    for packet in range(state.ngaussian):
        moments = _cross_moments_v252(
            q_c, p_c, width_c, chirp_c,
            state.q[packet], state.p[packet], state.widths[packet], state.chirps[packet],
            maximum_order=2,
        )
        kinetic = _kinetic_polynomial_v252(
            state.q[packet], state.p[packet], state.widths[packet],
            state.chirps[packet], model.mass_au,
        )
        coefficient = state.coefficients[packet]
        h_overlap += (
            _integrate_polynomial_v252(kinetic, moments) * coefficient
            + moments[0] * (model.H0 @ coefficient)
            + moments[1] * (model.H1 @ coefficient)
            + moments[2] * (model.H2 @ coefficient)
        )
    return dot_overlap + 1.0j * h_overlap


def _candidate_evaluation_data_v253(
    state, model, candidate, settings, *, metric_system=None
):
    state = state.validate()
    model = model.validate()
    settings = settings.validate()
    candidate = candidate.validate(packet_count=state.ngaussian)
    _validate_width_domain_v252(state, settings.tdvp_settings)
    if state.nstate != model.nstate:
        raise ValueError("spawn state/model electronic dimensions differ.")
    system = metric_system
    if system is None:
        system = build_adaptive_variational_metric_system_v252(
            state, model, settings=settings.tdvp_settings
        )
    else:
        if system.metric.shape != (state.parameter_count, state.parameter_count):
            raise ValueError("precomputed spawn metric has the wrong dimension.")
    geometry = (
        float(candidate.q), float(candidate.p),
        float(candidate.width), float(candidate.chirp),
    )
    raw_coupling = _residual_coupling_for_geometry_v253(
        state, model, system.velocity, geometry
    )
    current_couplings = np.zeros((state.ngaussian, state.nstate), dtype=complex)
    for index in range(state.ngaussian):
        current_couplings[index] = _residual_coupling_for_geometry_v253(
            state,
            model,
            system.velocity,
            (state.q[index], state.p[index], state.widths[index], state.chirps[index]),
        )
    candidate_to_current = _cross_overlap_v253(
        np.asarray([candidate.q]), np.asarray([candidate.p]),
        np.asarray([candidate.width]), np.asarray([candidate.chirp]), state,
    )[0]
    current_overlap = state.nuclear_overlap_matrix()
    current_bra_candidate = np.conj(candidate_to_current)
    projection_coefficients, _, _, _, _, _ = _svd_basis_solve_v253(
        current_overlap, current_bra_candidate[:, None], settings
    )
    projection_coefficients = projection_coefficients[:, 0]
    captured_norm = float(
        np.real(np.vdot(current_bra_candidate, projection_coefficients))
    )
    novelty = 1.0 - captured_norm
    if novelty < -settings.structural_tolerance:
        raise ValueError("spawn novelty became negative beyond structural tolerance.")
    novelty = max(float(novelty), 0.0)
    orthogonalized_coupling = raw_coupling - np.einsum(
        "i,ia->a", np.conj(projection_coefficients), current_couplings
    )
    residual_capture = float(
        np.linalg.norm(orthogonalized_coupling)
        / np.sqrt(max(novelty, np.finfo(float).tiny))
    )
    q = np.concatenate((state.q, [candidate.q]))
    p = np.concatenate((state.p, [candidate.p]))
    widths = np.concatenate((state.widths, [candidate.width]))
    chirps = np.concatenate((state.chirps, [candidate.chirp]))
    enlarged_overlap = _nuclear_overlap_from_arrays_v253(q, p, widths, chirps)
    enlarged_singular_values = np.linalg.svd(
        enlarged_overlap, compute_uv=False
    )
    largest = float(enlarged_singular_values[0])
    cutoff = max(
        settings.projection_absolute_cutoff,
        settings.projection_relative_cutoff * largest,
    )
    rank = int(np.count_nonzero(enlarged_singular_values > cutoff))
    condition = (
        float(enlarged_singular_values[0] / enlarged_singular_values[-1])
        if enlarged_singular_values[-1] > 0.0
        else float("inf")
    )
    reasons = []
    if state.ngaussian >= settings.maximum_packet_count:
        reasons.append("maximum-packet-count")
    if candidate.width < settings.tdvp_settings.minimum_width or candidate.width > settings.tdvp_settings.maximum_width:
        reasons.append("width-domain")
    if abs(candidate.chirp) > settings.tdvp_settings.maximum_absolute_chirp:
        reasons.append("chirp-domain")
    if novelty < settings.minimum_spawn_novelty:
        reasons.append("insufficient-novelty")
    if rank != state.ngaussian + 1:
        reasons.append("rank-deficient-enlarged-basis")
    if condition > settings.maximum_basis_condition_number:
        reasons.append("ill-conditioned-enlarged-basis")
    if residual_capture < settings.spawn_residual_capture_threshold:
        reasons.append("residual-score-below-threshold")
    return {
        "metric_velocity": system.velocity,
        "raw_coupling": raw_coupling,
        "current_couplings": current_couplings,
        "projection_coefficients": projection_coefficients,
        "orthogonalized_coupling": orthogonalized_coupling,
        "novelty": novelty,
        "maximum_existing_overlap": float(np.max(np.abs(candidate_to_current))),
        "enlarged_singular_values": enlarged_singular_values,
        "enlarged_cutoff": cutoff,
        "enlarged_rank": rank,
        "enlarged_condition": condition,
        "residual_capture": residual_capture,
        "rejection_reasons": tuple(reasons),
        "admitted": not reasons,
    }


@dataclass(frozen=True)
class SpawnCandidateEvaluationV253:
    state: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV252
    settings: ControlledBasisSettingsV253
    candidate: SpawnCandidateV253
    metric_velocity: np.ndarray
    raw_residual_coupling: np.ndarray
    current_basis_residual_couplings: np.ndarray
    projection_coefficients: np.ndarray
    orthogonalized_residual_coupling: np.ndarray
    novelty: float
    maximum_existing_overlap: float
    enlarged_singular_values: np.ndarray
    enlarged_cutoff: float
    enlarged_rank: int
    enlarged_condition_number: float
    residual_capture: float
    rejection_reasons: tuple
    admitted: bool

    def __post_init__(self):
        for name in (
            "metric_velocity", "raw_residual_coupling",
            "current_basis_residual_couplings", "projection_coefficients",
            "orthogonalized_residual_coupling",
        ):
            dtype = float if name == "metric_velocity" else complex
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=dtype).copy())
        object.__setattr__(
            self, "enlarged_singular_values",
            np.asarray(self.enlarged_singular_values, dtype=float).copy(),
        )
        object.__setattr__(self, "rejection_reasons", tuple(self.rejection_reasons))

    def validate(self):
        expected = _candidate_evaluation_data_v253(
            self.state, self.model, self.candidate, self.settings
        )
        tolerance = self.settings.structural_tolerance
        arrays = {
            "metric_velocity": "metric_velocity",
            "raw_residual_coupling": "raw_coupling",
            "current_basis_residual_couplings": "current_couplings",
            "projection_coefficients": "projection_coefficients",
            "orthogonalized_residual_coupling": "orthogonalized_coupling",
            "enlarged_singular_values": "enlarged_singular_values",
        }
        for stored_name, expected_name in arrays.items():
            if _scaled_norm_v252(getattr(self, stored_name), expected[expected_name]) > tolerance:
                raise ValueError(f"stored spawn evaluation {stored_name} is inconsistent.")
        scalars = {
            "novelty": "novelty",
            "maximum_existing_overlap": "maximum_existing_overlap",
            "enlarged_cutoff": "enlarged_cutoff",
            "enlarged_condition_number": "enlarged_condition",
            "residual_capture": "residual_capture",
        }
        for stored_name, expected_name in scalars.items():
            left = float(getattr(self, stored_name))
            right = float(expected[expected_name])
            if not (np.isinf(left) and np.isinf(right)) and abs(left - right) > tolerance:
                raise ValueError(f"stored spawn evaluation {stored_name} is inconsistent.")
        if int(self.enlarged_rank) != expected["enlarged_rank"]:
            raise ValueError("stored spawn enlarged rank is inconsistent.")
        if tuple(self.rejection_reasons) != expected["rejection_reasons"]:
            raise ValueError("stored spawn rejection reasons are inconsistent.")
        if type(self.admitted) is not bool or self.admitted != expected["admitted"]:
            raise ValueError("stored spawn admission decision is inconsistent.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "candidate": self.candidate.as_dict(),
            "state_fingerprint": _sha256_v252(self.state.as_dict()),
            "model_fingerprint": self.model.fingerprint(),
            "metric_velocity": self.metric_velocity.tolist(),
            "raw_residual_coupling": _complex_pairs_v252(self.raw_residual_coupling),
            "current_basis_residual_couplings": _complex_pairs_v252(
                self.current_basis_residual_couplings
            ),
            "projection_coefficients": _complex_pairs_v252(self.projection_coefficients),
            "orthogonalized_residual_coupling": _complex_pairs_v252(
                self.orthogonalized_residual_coupling
            ),
            "novelty": float(self.novelty),
            "maximum_existing_overlap": float(self.maximum_existing_overlap),
            "enlarged_singular_values": self.enlarged_singular_values.tolist(),
            "enlarged_cutoff": float(self.enlarged_cutoff),
            "enlarged_rank": int(self.enlarged_rank),
            "enlarged_condition_number": float(self.enlarged_condition_number),
            "residual_capture": float(self.residual_capture),
            "rejection_reasons": list(self.rejection_reasons),
            "admitted": self.admitted,
        }


def evaluate_spawn_candidate_v253(
    state, model, candidate, *, settings=ControlledBasisSettingsV253(),
    _metric_system=None,
):
    data = _candidate_evaluation_data_v253(
        state, model, candidate, settings, metric_system=_metric_system
    )
    return SpawnCandidateEvaluationV253(
        state=state,
        model=model,
        settings=settings,
        candidate=candidate,
        metric_velocity=data["metric_velocity"],
        raw_residual_coupling=data["raw_coupling"],
        current_basis_residual_couplings=data["current_couplings"],
        projection_coefficients=data["projection_coefficients"],
        orthogonalized_residual_coupling=data["orthogonalized_coupling"],
        novelty=data["novelty"],
        maximum_existing_overlap=data["maximum_existing_overlap"],
        enlarged_singular_values=data["enlarged_singular_values"],
        enlarged_cutoff=data["enlarged_cutoff"],
        enlarged_rank=data["enlarged_rank"],
        enlarged_condition_number=data["enlarged_condition"],
        residual_capture=data["residual_capture"],
        rejection_reasons=data["rejection_reasons"],
        admitted=data["admitted"],
    ).validate()


def _projection_after_removal_v253(state, model, remove_index, event_kind, settings):
    keep = np.asarray([index for index in range(state.ngaussian) if index != remove_index])
    return project_adaptive_state_v253(
        state,
        model,
        q=state.q[keep],
        p=state.p[keep],
        widths=state.widths[keep],
        chirps=state.chirps[keep],
        event_kind=event_kind,
        settings=settings,
    )


@dataclass(frozen=True)
class BasisLifecycleEventV253:
    event_kind: str
    before: ThawedGaussianSpinorStateV252
    after: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV252
    settings: ControlledBasisSettingsV253
    packet_ids_before: tuple
    packet_ids_after: tuple
    packet_ages_before: tuple
    packet_ages_after: tuple
    next_packet_serial_before: int
    next_packet_serial_after: int
    candidate_evaluations: tuple = field(default_factory=tuple)
    selected_candidate: SpawnCandidateV253 | None = None
    projection: BasisProjectionReceiptV253 | None = None
    added_packet_id: str | None = None
    removed_packet_id: str | None = None
    reason: str = ""

    def validate(self):
        self.before.validate()
        self.after.validate()
        self.model.validate()
        self.settings.validate()
        ids_before, ages_before = _validate_packet_metadata_v253(
            self.packet_ids_before, self.packet_ages_before, self.before.ngaussian
        )
        ids_after, ages_after = _validate_packet_metadata_v253(
            self.packet_ids_after, self.packet_ages_after, self.after.ngaussian
        )
        serial_before = _positive_integer_v253(
            self.next_packet_serial_before, "next_packet_serial_before", allow_zero=True
        )
        serial_after = _positive_integer_v253(
            self.next_packet_serial_after, "next_packet_serial_after", allow_zero=True
        )
        if self.event_kind not in ("none", "spawn", "prune", "merge"):
            raise ValueError("basis lifecycle event kind is invalid.")
        if not isinstance(self.reason, str) or not self.reason:
            raise ValueError("basis lifecycle event reason is missing.")
        for evaluation in self.candidate_evaluations:
            evaluation.validate()
            if _state_distance_v253(evaluation.state, self.before) > self.settings.structural_tolerance:
                raise ValueError("spawn evaluation is not bound to the event input.")
        if self.event_kind == "none":
            if _state_distance_v253(self.before, self.after) > self.settings.structural_tolerance:
                raise ValueError("no-op lifecycle event changed the state.")
            if ids_before != ids_after or ages_before != ages_after:
                raise ValueError("no-op lifecycle event changed packet metadata.")
            if any(value is not None for value in (
                self.selected_candidate, self.projection,
                self.added_packet_id, self.removed_packet_id,
            )):
                raise ValueError("no-op lifecycle event contains topology data.")
            if serial_after != serial_before:
                raise ValueError("no-op lifecycle event changed the packet serial.")
            return self
        if self.projection is None:
            raise ValueError("topology event is missing its projection receipt.")
        self.projection.validate()
        if self.projection.event_kind != self.event_kind:
            raise ValueError("topology event/projection kinds disagree.")
        if _state_distance_v253(self.projection.source, self.before) > self.settings.structural_tolerance:
            raise ValueError("topology projection source is inconsistent.")
        if _state_distance_v253(self.projection.projected, self.after) > self.settings.structural_tolerance:
            raise ValueError("topology projection target is inconsistent.")
        if self.event_kind == "spawn":
            if self.selected_candidate is None or self.added_packet_id is None or self.removed_packet_id is not None:
                raise ValueError("spawn event metadata is incomplete.")
            self.selected_candidate.validate(packet_count=self.before.ngaussian)
            if len(ids_after) != len(ids_before) + 1 or ids_after[:-1] != ids_before:
                raise ValueError("spawn packet ID transition is inconsistent.")
            expected_id = f"g{serial_before:06d}"
            if self.added_packet_id != expected_id or ids_after[-1] != expected_id:
                raise ValueError("spawn packet ID does not match the stable serial.")
            if ages_after[:-1] != ages_before or ages_after[-1] != 0:
                raise ValueError("spawn packet ages are inconsistent.")
            if serial_after != serial_before + 1:
                raise ValueError("spawn serial did not advance exactly once.")
            admitted = [item for item in self.candidate_evaluations if item.admitted]
            if not admitted:
                raise ValueError("spawn event has no admitted candidate.")
            expected_selection = sorted(
                admitted,
                key=lambda item: (-item.residual_capture, item.candidate.canonical_key()),
            )[0].candidate
            if self.selected_candidate.as_dict() != expected_selection.as_dict():
                raise ValueError("spawn event did not select the highest admitted score.")
        else:
            if self.selected_candidate is not None or self.added_packet_id is not None:
                raise ValueError("removal event contains spawn metadata.")
            if self.removed_packet_id not in ids_before:
                raise ValueError("removed packet ID was absent before the event.")
            remove_index = ids_before.index(self.removed_packet_id)
            expected_ids = ids_before[:remove_index] + ids_before[remove_index + 1 :]
            expected_ages = ages_before[:remove_index] + ages_before[remove_index + 1 :]
            if ids_after != expected_ids or ages_after != expected_ages:
                raise ValueError("removal event packet metadata is inconsistent.")
            if serial_after != serial_before:
                raise ValueError("removal event changed the packet serial.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": CONTROLLED_BASIS_EVENT_SCHEMA_V253,
            "event_kind": self.event_kind,
            "before": self.before.as_dict(),
            "after": self.after.as_dict(),
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "packet_ids_before": list(self.packet_ids_before),
            "packet_ids_after": list(self.packet_ids_after),
            "packet_ages_before": list(self.packet_ages_before),
            "packet_ages_after": list(self.packet_ages_after),
            "next_packet_serial_before": int(self.next_packet_serial_before),
            "next_packet_serial_after": int(self.next_packet_serial_after),
            "candidate_evaluations": [item.as_dict() for item in self.candidate_evaluations],
            "selected_candidate": (
                None if self.selected_candidate is None else self.selected_candidate.as_dict()
            ),
            "projection": None if self.projection is None else self.projection.as_dict(),
            "added_packet_id": self.added_packet_id,
            "removed_packet_id": self.removed_packet_id,
            "reason": self.reason,
        }

    def fingerprint(self):
        return _sha256_v252(self.as_dict())


def _no_event_v253(
    state, model, settings, packet_ids, packet_ages, next_packet_serial,
    *, evaluations=(), reason,
):
    return BasisLifecycleEventV253(
        event_kind="none",
        before=state,
        after=state,
        model=model,
        settings=settings,
        packet_ids_before=tuple(packet_ids),
        packet_ids_after=tuple(packet_ids),
        packet_ages_before=tuple(packet_ages),
        packet_ages_after=tuple(packet_ages),
        next_packet_serial_before=next_packet_serial,
        next_packet_serial_after=next_packet_serial,
        candidate_evaluations=tuple(evaluations),
        reason=reason,
    ).validate()


def adapt_basis_once_v253(
    state,
    model,
    *,
    packet_ids,
    packet_ages,
    next_packet_serial,
    settings=ControlledBasisSettingsV253(),
    spawn_candidates=None,
):
    """Apply the frozen merge -> prune -> spawn policy at one checkpoint."""

    state = state.validate()
    model = model.validate()
    settings = settings.validate()
    _validate_width_domain_v252(state, settings.tdvp_settings)
    packet_ids, packet_ages = _validate_packet_metadata_v253(
        packet_ids, packet_ages, state.ngaussian
    )
    next_packet_serial = _positive_integer_v253(
        next_packet_serial, "next_packet_serial", allow_zero=True
    )
    if not settings.minimum_packet_count <= state.ngaussian <= settings.maximum_packet_count:
        raise ValueError("current Gaussian count is outside the lifecycle domain.")
    if state.nstate != model.nstate:
        raise ValueError("lifecycle state/model electronic dimensions differ.")

    if state.ngaussian > settings.minimum_packet_count:
        overlap = state.nuclear_overlap_matrix()
        merge_options = []
        for left in range(state.ngaussian):
            for right in range(left + 1, state.ngaussian):
                if min(packet_ages[left], packet_ages[right]) < settings.minimum_packet_age_steps:
                    continue
                pair_overlap = float(abs(overlap[left, right]))
                if pair_overlap < settings.minimum_merge_overlap:
                    continue
                for remove_index in (left, right):
                    try:
                        projection = _projection_after_removal_v253(
                            state, model, remove_index, "merge", settings
                        )
                    except ValueError:
                        continue
                    merge_options.append(
                        (
                            projection.relative_projection_loss,
                            abs(projection.energy_jump_hartree),
                            packet_ids[remove_index],
                            -pair_overlap,
                            remove_index,
                            projection,
                        )
                    )
        if merge_options:
            _, _, _, _, remove_index, projection = sorted(merge_options)[0]
            removed = packet_ids[remove_index]
            ids_after = packet_ids[:remove_index] + packet_ids[remove_index + 1 :]
            ages_after = packet_ages[:remove_index] + packet_ages[remove_index + 1 :]
            return BasisLifecycleEventV253(
                event_kind="merge",
                before=state,
                after=projection.projected,
                model=model,
                settings=settings,
                packet_ids_before=packet_ids,
                packet_ids_after=ids_after,
                packet_ages_before=packet_ages,
                packet_ages_after=ages_after,
                next_packet_serial_before=next_packet_serial,
                next_packet_serial_after=next_packet_serial,
                projection=projection,
                removed_packet_id=removed,
                reason="accepted highest-overlap pair with minimum projection loss",
            ).validate()

        prune_options = []
        populations = np.sum(np.abs(state.coefficients) ** 2, axis=1)
        for remove_index in range(state.ngaussian):
            if packet_ages[remove_index] < settings.minimum_packet_age_steps:
                continue
            activation_timed_out = (
                populations[remove_index] < settings.shape_activation_population
                and packet_ages[remove_index] >= settings.maximum_activation_age_steps
            )
            if (
                populations[remove_index] > settings.maximum_prune_coefficient_population
                and not activation_timed_out
            ):
                continue
            if (
                populations[remove_index] < settings.shape_activation_population
                and packet_ages[remove_index] < settings.maximum_activation_age_steps
            ):
                continue
            try:
                projection = _projection_after_removal_v253(
                    state, model, remove_index, "prune", settings
                )
            except ValueError:
                continue
            prune_options.append(
                (
                    projection.relative_projection_loss,
                    abs(projection.energy_jump_hartree),
                    float(populations[remove_index]),
                    packet_ids[remove_index],
                    remove_index,
                    projection,
                )
            )
        if prune_options:
            _, _, _, _, remove_index, projection = sorted(prune_options)[0]
            removed = packet_ids[remove_index]
            ids_after = packet_ids[:remove_index] + packet_ids[remove_index + 1 :]
            ages_after = packet_ages[:remove_index] + packet_ages[remove_index + 1 :]
            return BasisLifecycleEventV253(
                event_kind="prune",
                before=state,
                after=projection.projected,
                model=model,
                settings=settings,
                packet_ids_before=packet_ids,
                packet_ids_after=ids_after,
                packet_ages_before=packet_ages,
                packet_ages_after=ages_after,
                next_packet_serial_before=next_packet_serial,
                next_packet_serial_after=next_packet_serial,
                projection=projection,
                removed_packet_id=removed,
                reason="accepted lowest-loss age-eligible low-population removal",
            ).validate()

    if state.ngaussian >= settings.maximum_packet_count:
        return _no_event_v253(
            state, model, settings, packet_ids, packet_ages, next_packet_serial,
            reason="maximum packet count reached",
        )
    populations = np.sum(np.abs(state.coefficients) ** 2, axis=1)
    if np.any(populations < settings.shape_activation_population):
        return _no_event_v253(
            state, model, settings, packet_ids, packet_ages, next_packet_serial,
            reason="dormant packet coefficient activation is still in progress",
        )
    if spawn_candidates is None:
        spawn_candidates = generate_spawn_candidates_v253(state, settings=settings)
    else:
        spawn_candidates = tuple(
            candidate.validate(packet_count=state.ngaussian)
            for candidate in spawn_candidates
        )
    spawn_metric_system = build_adaptive_variational_metric_system_v252(
        state, model, settings=settings.tdvp_settings
    )
    evaluations = tuple(
        evaluate_spawn_candidate_v253(
            state, model, candidate, settings=settings,
            _metric_system=spawn_metric_system,
        )
        for candidate in spawn_candidates
    )
    admitted = [evaluation for evaluation in evaluations if evaluation.admitted]
    if not admitted:
        return _no_event_v253(
            state, model, settings, packet_ids, packet_ages, next_packet_serial,
            evaluations=evaluations,
            reason="no residual candidate passed every spawn gate",
        )
    selected = sorted(
        admitted,
        key=lambda item: (-item.residual_capture, item.candidate.canonical_key()),
    )[0]
    candidate = selected.candidate
    q = np.concatenate((state.q, [candidate.q]))
    p = np.concatenate((state.p, [candidate.p]))
    widths = np.concatenate((state.widths, [candidate.width]))
    chirps = np.concatenate((state.chirps, [candidate.chirp]))
    projection = project_adaptive_state_v253(
        state,
        model,
        q=q,
        p=p,
        widths=widths,
        chirps=chirps,
        event_kind="spawn",
        settings=settings,
    )
    added = f"g{next_packet_serial:06d}"
    if added in packet_ids:
        raise ValueError("next packet serial collides with an existing stable ID.")
    return BasisLifecycleEventV253(
        event_kind="spawn",
        before=state,
        after=projection.projected,
        model=model,
        settings=settings,
        packet_ids_before=packet_ids,
        packet_ids_after=packet_ids + (added,),
        packet_ages_before=packet_ages,
        packet_ages_after=packet_ages + (0,),
        next_packet_serial_before=next_packet_serial,
        next_packet_serial_after=next_packet_serial + 1,
        candidate_evaluations=evaluations,
        selected_candidate=candidate,
        projection=projection,
        added_packet_id=added,
        reason="accepted highest residual-capture candidate after novelty and conditioning gates",
    ).validate()


def _raw_adaptive_metric_rhs_v253(state, model, tdvp_settings):
    """Build the exact v0.25.2 metric/RHS without prematurely solving null shapes."""

    state = _validate_width_domain_v252(state, tdvp_settings)
    model = model.validate()
    if state.nstate != model.nstate:
        raise ValueError("controlled metric state/model dimensions differ.")
    terms = _tangent_terms_v252(state)
    count = len(terms)
    metric = np.zeros((count, count), dtype=float)
    rhs = np.zeros(count, dtype=float)
    moments = {}
    for i in range(state.ngaussian):
        for j in range(state.ngaussian):
            moments[(i, j)] = _cross_moments_v252(
                state.q[i], state.p[i], state.widths[i], state.chirps[i],
                state.q[j], state.p[j], state.widths[j], state.chirps[j],
                maximum_order=4,
            )
    identity = np.eye(model.nstate, dtype=complex)
    for mu, (i, vector_i, polynomial_i) in enumerate(terms):
        conjugate_i = np.conj(polynomial_i)
        for nu, (j, vector_j, polynomial_j) in enumerate(terms):
            polynomial = np.convolve(conjugate_i, polynomial_j)
            nuclear = _integrate_polynomial_v252(polynomial, moments[(i, j)])
            metric[mu, nu] = float(np.real(np.vdot(vector_i, vector_j) * nuclear))
        projection = 0.0 + 0.0j
        for j in range(state.ngaussian):
            kinetic = _kinetic_polynomial_v252(
                state.q[j], state.p[j], state.widths[j], state.chirps[j], model.mass_au
            )
            operators = (
                kinetic[0] * identity + model.H0,
                kinetic[1] * identity + model.H1,
                kinetic[2] * identity + model.H2,
            )
            for bra_degree, bra_value in enumerate(conjugate_i):
                for operator_degree, operator in enumerate(operators):
                    projection += (
                        bra_value
                        * moments[(i, j)][bra_degree + operator_degree]
                        * np.vdot(vector_i, operator @ state.coefficients[j])
                    )
        rhs[mu] = float(np.imag(projection))
    return 0.5 * (metric + metric.T), rhs


def _active_parameter_indices_v253(state, active_shape_mask):
    mask = np.asarray(active_shape_mask)
    if mask.shape != (state.ngaussian,) or mask.dtype != bool:
        raise ValueError("active shape mask must be a Boolean per-packet vector.")
    block = state.ngaussian * state.nstate
    indices = list(range(2 * block))
    for offset in range(4):
        start = 2 * block + offset * state.ngaussian
        indices.extend(start + index for index in range(state.ngaussian) if mask[index])
    return np.asarray(indices, dtype=int)


@dataclass(frozen=True)
class ControlledMetricSystemV253:
    state: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV252
    tdvp_settings: AdaptiveVariationalSettingsV252
    active_shape_mask: np.ndarray
    active_parameter_indices: np.ndarray
    full_metric: np.ndarray
    full_rhs: np.ndarray
    reduced_metric: np.ndarray
    reduced_rhs: np.ndarray
    reduced_velocity: np.ndarray
    full_velocity: np.ndarray
    solve_receipt: MetricSolveReceiptV251

    def __post_init__(self):
        object.__setattr__(self, "active_shape_mask", np.asarray(self.active_shape_mask, dtype=bool).copy())
        object.__setattr__(self, "active_parameter_indices", np.asarray(self.active_parameter_indices, dtype=int).copy())
        for name in (
            "full_metric", "full_rhs", "reduced_metric", "reduced_rhs",
            "reduced_velocity", "full_velocity",
        ):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())

    def validate(self):
        self.state.validate(require_normalized=False)
        self.model.validate()
        self.tdvp_settings.validate()
        expected_indices = _active_parameter_indices_v253(self.state, self.active_shape_mask)
        if not np.array_equal(self.active_parameter_indices, expected_indices):
            raise ValueError("controlled metric active parameter indices are inconsistent.")
        metric, rhs = _raw_adaptive_metric_rhs_v253(
            self.state, self.model, self.tdvp_settings
        )
        indices = expected_indices
        reduced_metric = metric[np.ix_(indices, indices)]
        reduced_rhs = rhs[indices]
        velocity, receipt = solve_variational_metric_v251(
            reduced_metric, reduced_rhs, settings=self.tdvp_settings
        )
        full_velocity = np.zeros(self.state.parameter_count, dtype=float)
        full_velocity[indices] = velocity
        tolerance = self.tdvp_settings.structural_tolerance
        expected_arrays = {
            "full_metric": metric,
            "full_rhs": rhs,
            "reduced_metric": reduced_metric,
            "reduced_rhs": reduced_rhs,
            "reduced_velocity": velocity,
            "full_velocity": full_velocity,
        }
        for name, expected in expected_arrays.items():
            if _scaled_norm_v252(getattr(self, name), expected) > tolerance:
                raise ValueError(f"stored controlled metric {name} is inconsistent.")
        if _scaled_norm_v252(self.solve_receipt.singular_values, receipt.singular_values) > tolerance:
            raise ValueError("stored controlled metric spectrum is inconsistent.")
        if self.solve_receipt.rank != receipt.rank or self.solve_receipt.nullity != receipt.nullity:
            raise ValueError("stored controlled metric rank/nullity is inconsistent.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "state_fingerprint": _sha256_v252(self.state.as_dict()),
            "model_fingerprint": self.model.fingerprint(),
            "active_shape_mask": self.active_shape_mask.tolist(),
            "active_parameter_indices": self.active_parameter_indices.tolist(),
            "full_metric": self.full_metric.tolist(),
            "full_rhs": self.full_rhs.tolist(),
            "reduced_metric": self.reduced_metric.tolist(),
            "reduced_rhs": self.reduced_rhs.tolist(),
            "reduced_velocity": self.reduced_velocity.tolist(),
            "full_velocity": self.full_velocity.tolist(),
            "solve_receipt": self.solve_receipt.as_dict(),
        }


def build_controlled_metric_system_v253(
    state, model, active_shape_mask, *, tdvp_settings=AdaptiveVariationalSettingsV252()
):
    state = state.validate(require_normalized=False)
    model = model.validate()
    tdvp_settings = tdvp_settings.validate()
    mask = np.asarray(active_shape_mask, dtype=bool)
    indices = _active_parameter_indices_v253(state, mask)
    metric, rhs = _raw_adaptive_metric_rhs_v253(state, model, tdvp_settings)
    reduced_metric = metric[np.ix_(indices, indices)]
    reduced_rhs = rhs[indices]
    velocity, receipt = solve_variational_metric_v251(
        reduced_metric, reduced_rhs, settings=tdvp_settings
    )
    full_velocity = np.zeros(state.parameter_count, dtype=float)
    full_velocity[indices] = velocity
    return ControlledMetricSystemV253(
        state=state,
        model=model,
        tdvp_settings=tdvp_settings,
        active_shape_mask=mask,
        active_parameter_indices=indices,
        full_metric=metric,
        full_rhs=rhs,
        reduced_metric=reduced_metric,
        reduced_rhs=reduced_rhs,
        reduced_velocity=velocity,
        full_velocity=full_velocity,
        solve_receipt=receipt,
    ).validate()


@dataclass(frozen=True)
class CoefficientActivationStepV253:
    start: ThawedGaussianSpinorStateV252
    end: ThawedGaussianSpinorStateV252
    dt_au: float
    model: QuadraticSpinHamiltonianV252
    settings: AdaptiveVariationalSettingsV252
    active_shape_mask: np.ndarray
    active_parameter_indices: np.ndarray
    midpoint_parameters: np.ndarray
    midpoint_system: ControlledMetricSystemV253
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
        object.__setattr__(self, "active_shape_mask", np.asarray(self.active_shape_mask, dtype=bool).copy())
        object.__setattr__(self, "active_parameter_indices", np.asarray(self.active_parameter_indices, dtype=int).copy())
        for name in ("midpoint_parameters", "nonlinear_residual"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())

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
        self.model.validate()
        tolerance = max(settings.structural_tolerance, settings.maximum_step_norm_drift)
        self.start.validate(tolerance=tolerance)
        self.end.validate(tolerance=tolerance)
        if self.start.ngaussian != self.end.ngaussian or self.start.nstate != self.end.nstate:
            raise ValueError("coefficient-activation step changed basis dimension.")
        dt = float(self.dt_au)
        if not np.isfinite(dt) or dt == 0.0:
            raise ValueError("coefficient-activation dt must be finite and nonzero.")
        if abs(self.end.time_au - self.start.time_au - dt) > settings.structural_tolerance:
            raise ValueError("coefficient-activation endpoint time is inconsistent.")
        indices = _active_parameter_indices_v253(self.start, self.active_shape_mask)
        if not np.array_equal(indices, self.active_parameter_indices):
            raise ValueError("coefficient-activation indices are inconsistent.")
        theta_start = pack_adaptive_variational_parameters_v252(self.start)
        theta_end = pack_adaptive_variational_parameters_v252(self.end)
        frozen = np.ones(self.start.parameter_count, dtype=bool)
        frozen[indices] = False
        if np.max(np.abs(theta_end[frozen] - theta_start[frozen]), initial=0.0) > settings.structural_tolerance:
            raise ValueError("inactive newborn shape coordinates moved.")
        expected_midpoint = 0.5 * (theta_start + theta_end)
        if _scaled_norm_v252(self.midpoint_parameters, expected_midpoint) > settings.structural_tolerance:
            raise ValueError("coefficient-activation midpoint parameters are inconsistent.")
        midpoint_state = state_from_adaptive_variational_parameters_v252(
            expected_midpoint,
            ngaussian=self.start.ngaussian,
            nstate=self.start.nstate,
            time_au=self.start.time_au + 0.5 * dt,
        )
        expected_system = build_controlled_metric_system_v253(
            midpoint_state, self.model, self.active_shape_mask, tdvp_settings=settings
        )
        self.midpoint_system.validate()
        for name in ("reduced_metric", "reduced_rhs", "reduced_velocity", "full_velocity"):
            if _scaled_norm_v252(getattr(self.midpoint_system, name), getattr(expected_system, name)) > settings.structural_tolerance:
                raise ValueError(f"stored activation midpoint {name} is inconsistent.")
        expected_residual = theta_end[indices] - theta_start[indices] - dt * expected_system.reduced_velocity
        if _scaled_norm_v252(self.nonlinear_residual, expected_residual) > settings.structural_tolerance:
            raise ValueError("stored activation nonlinear residual is inconsistent.")
        expected_norm = float(np.linalg.norm(expected_residual))
        if abs(self.nonlinear_residual_norm - expected_norm) > settings.structural_tolerance:
            raise ValueError("stored activation nonlinear residual norm is inconsistent.")
        if expected_norm > settings.nonlinear_residual_tolerance:
            raise ValueError("activation nonlinear residual exceeds tolerance.")
        if type(self.nonlinear_success) is not bool or not self.nonlinear_success or int(self.nonlinear_status) <= 0:
            raise ValueError("activation nonlinear solver did not report success.")
        expected_scalars = {
            "start_norm": self.start.generalized_norm,
            "end_norm": self.end.generalized_norm,
            "start_energy_hartree": adaptive_variational_energy_v252(self.start, self.model),
            "end_energy_hartree": adaptive_variational_energy_v252(self.end, self.model),
        }
        for name, expected in expected_scalars.items():
            if abs(float(getattr(self, name)) - float(expected)) > settings.structural_tolerance:
                raise ValueError(f"stored activation {name} is inconsistent.")
        if abs(self.norm_change) > settings.maximum_step_norm_drift:
            raise ValueError("activation step norm drift exceeds the release gate.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "mode": "coefficient-only activation for dormant packet shapes",
            "start": self.start.as_dict(),
            "end": self.end.as_dict(),
            "dt_au": float(self.dt_au),
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "active_shape_mask": self.active_shape_mask.tolist(),
            "active_parameter_indices": self.active_parameter_indices.tolist(),
            "midpoint_parameters": self.midpoint_parameters.tolist(),
            "midpoint_system": self.midpoint_system.as_dict(),
            "nonlinear": {
                "success": self.nonlinear_success,
                "status": int(self.nonlinear_status),
                "message": self.nonlinear_message,
                "function_evaluations": int(self.nonlinear_function_evaluations),
                "jacobian_evaluations": self.nonlinear_jacobian_evaluations,
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


def coefficient_activation_implicit_step_v253(
    state, model, dt_au, active_shape_mask, *, settings=AdaptiveVariationalSettingsV252()
):
    settings = settings.validate()
    model = model.validate()
    state = state.validate(
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift)
    )
    mask = np.asarray(active_shape_mask, dtype=bool)
    indices = _active_parameter_indices_v253(state, mask)
    if np.all(mask):
        raise ValueError("coefficient-activation step requires at least one dormant shape.")
    dt = float(dt_au)
    if not np.isfinite(dt) or dt == 0.0:
        raise ValueError("coefficient-activation dt must be finite and nonzero.")
    theta_start = pack_adaptive_variational_parameters_v252(state)
    initial_system = build_controlled_metric_system_v253(
        state, model, mask, tdvp_settings=settings
    )
    predictor = theta_start[indices] + dt * initial_system.reduced_velocity

    def residual(active_endpoint):
        theta_end = theta_start.copy()
        theta_end[indices] = np.asarray(active_endpoint, dtype=float)
        midpoint = 0.5 * (theta_start + theta_end)
        midpoint_state = state_from_adaptive_variational_parameters_v252(
            midpoint,
            ngaussian=state.ngaussian,
            nstate=state.nstate,
            time_au=state.time_au + 0.5 * dt,
        )
        system = build_controlled_metric_system_v253(
            midpoint_state, model, mask, tdvp_settings=settings
        )
        return np.asarray(active_endpoint) - theta_start[indices] - dt * system.reduced_velocity

    predictor_residual_norm = float(np.linalg.norm(residual(predictor)))
    solution = root(
        residual,
        predictor,
        method="hybr",
        options={"xtol": settings.nonlinear_xtol, "maxfev": settings.nonlinear_max_function_evaluations},
    )
    final_residual = np.asarray(residual(solution.x), dtype=float)
    final_residual_norm = float(np.linalg.norm(final_residual))
    if not bool(solution.success) or final_residual_norm > settings.nonlinear_residual_tolerance:
        raise RuntimeError(
            "coefficient-activation implicit midpoint solve failed: "
            f"success={bool(solution.success)}, status={int(solution.status)}, "
            f"residual={final_residual_norm:.6e}, message={solution.message}"
        )
    theta_end = theta_start.copy()
    theta_end[indices] = solution.x
    end = state_from_adaptive_variational_parameters_v252(
        theta_end,
        ngaussian=state.ngaussian,
        nstate=state.nstate,
        time_au=state.time_au + dt,
    ).validate(
        tolerance=max(settings.structural_tolerance, settings.maximum_step_norm_drift)
    )
    midpoint = 0.5 * (theta_start + theta_end)
    midpoint_state = state_from_adaptive_variational_parameters_v252(
        midpoint,
        ngaussian=state.ngaussian,
        nstate=state.nstate,
        time_au=state.time_au + 0.5 * dt,
    )
    midpoint_system = build_controlled_metric_system_v253(
        midpoint_state, model, mask, tdvp_settings=settings
    )
    return CoefficientActivationStepV253(
        start=state,
        end=end,
        dt_au=dt,
        model=model,
        settings=settings,
        active_shape_mask=mask,
        active_parameter_indices=indices,
        midpoint_parameters=midpoint,
        midpoint_system=midpoint_system,
        nonlinear_success=bool(solution.success),
        nonlinear_status=int(solution.status),
        nonlinear_message=str(solution.message),
        nonlinear_function_evaluations=int(solution.nfev),
        nonlinear_jacobian_evaluations=None if not hasattr(solution, "njev") else int(solution.njev),
        nonlinear_residual=final_residual,
        nonlinear_residual_norm=final_residual_norm,
        predictor_residual_norm=predictor_residual_norm,
        start_norm=state.generalized_norm,
        end_norm=end.generalized_norm,
        start_energy_hartree=adaptive_variational_energy_v252(state, model),
        end_energy_hartree=adaptive_variational_energy_v252(end, model),
    ).validate()


def controlled_tdvp_step_v253(state, model, dt_au, *, settings=ControlledBasisSettingsV253()):
    """Use v0.25.2 directly unless a zero/small-amplitude packet needs activation."""

    settings = settings.validate()
    state = state.validate()
    populations = np.sum(np.abs(state.coefficients) ** 2, axis=1)
    active_shape_mask = populations >= settings.shape_activation_population
    if np.all(active_shape_mask):
        return adaptive_implicit_midpoint_tdvp_step_v252(
            state, model, dt_au, settings=settings.tdvp_settings
        )
    return coefficient_activation_implicit_step_v253(
        state, model, dt_au, active_shape_mask, settings=settings.tdvp_settings
    )


V253_CONTROLLED_BASIS_CLAIMS = {
    "adaptive_width_multigaussian_tdvp_inherited": True,
    "analytic_residual_candidate_scoring_validated": True,
    "controlled_residual_driven_spawning_validated": True,
    "full_svd_enlarged_basis_projection_validated": True,
    "projection_guarded_pruning_validated": True,
    "overlap_projection_guarded_merging_validated": True,
    "stable_packet_identity_and_age_gates_validated": True,
    "one_topology_event_per_checkpoint_validated": True,
    "coefficient_only_newborn_activation_validated": True,
    "constant_electronic_gauge_covariance_validated": True,
    "gaussian_permutation_covariance_validated": True,
    "general_aims_branching_validated": False,
    "multidimensional_spawning_validated": False,
    "full_correlated_width_matrices_validated": False,
    "coordinate_dependent_electronic_gauge_covariance_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
    "general_ab_initio_soc_dynamics_accuracy_validated": False,
}


@dataclass(frozen=True)
class ControlledBasisStepV253:
    step_index: int
    start: ThawedGaussianSpinorStateV252
    end: ThawedGaussianSpinorStateV252
    packet_ids_start: tuple
    packet_ids_end: tuple
    packet_ages_start: tuple
    packet_ages_end: tuple
    next_packet_serial_start: int
    next_packet_serial_end: int
    tdvp_step: object
    lifecycle_event: BasisLifecycleEventV253
    settings: ControlledBasisSettingsV253

    def validate(self):
        index = _positive_integer_v253(self.step_index, "controlled basis step index")
        ids_start, ages_start = _validate_packet_metadata_v253(
            self.packet_ids_start, self.packet_ages_start, self.start.ngaussian
        )
        ids_end, ages_end = _validate_packet_metadata_v253(
            self.packet_ids_end, self.packet_ages_end, self.end.ngaussian
        )
        serial_start = _positive_integer_v253(
            self.next_packet_serial_start, "next_packet_serial_start", allow_zero=True
        )
        serial_end = _positive_integer_v253(
            self.next_packet_serial_end, "next_packet_serial_end", allow_zero=True
        )
        self.settings.validate()
        self.tdvp_step.validate()
        self.lifecycle_event.validate()
        if _state_distance_v253(self.start, self.tdvp_step.start) > self.settings.structural_tolerance:
            raise ValueError("controlled step is not bound to its TDVP start.")
        if _state_distance_v253(self.tdvp_step.end, self.lifecycle_event.before) > self.settings.structural_tolerance:
            raise ValueError("controlled step TDVP/event boundary is discontinuous.")
        if _state_distance_v253(self.lifecycle_event.after, self.end) > self.settings.structural_tolerance:
            raise ValueError("controlled step is not bound to its lifecycle endpoint.")
        if self.tdvp_step.settings.as_dict() != self.settings.tdvp_settings.as_dict():
            raise ValueError("controlled step changed the inherited TDVP settings.")
        if tuple(self.lifecycle_event.packet_ids_before) != ids_start:
            raise ValueError("controlled step changed IDs during fixed-basis propagation.")
        expected_event_ages = tuple(value + 1 for value in ages_start)
        if tuple(self.lifecycle_event.packet_ages_before) != expected_event_ages:
            raise ValueError("controlled step packet ages did not advance exactly once.")
        if tuple(self.lifecycle_event.packet_ids_after) != ids_end or tuple(self.lifecycle_event.packet_ages_after) != ages_end:
            raise ValueError("controlled step endpoint metadata is inconsistent.")
        if self.lifecycle_event.next_packet_serial_before != serial_start or self.lifecycle_event.next_packet_serial_after != serial_end:
            raise ValueError("controlled step packet serial transition is inconsistent.")
        if index % self.settings.adaptation_interval_steps != 0:
            if self.lifecycle_event.event_kind != "none" or self.lifecycle_event.reason != "adaptation interval not reached":
                raise ValueError("off-interval controlled step performed adaptation.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "step_index": int(self.step_index),
            "start": self.start.as_dict(),
            "end": self.end.as_dict(),
            "packet_ids_start": list(self.packet_ids_start),
            "packet_ids_end": list(self.packet_ids_end),
            "packet_ages_start": list(self.packet_ages_start),
            "packet_ages_end": list(self.packet_ages_end),
            "next_packet_serial_start": int(self.next_packet_serial_start),
            "next_packet_serial_end": int(self.next_packet_serial_end),
            "tdvp_step": self.tdvp_step.as_dict(),
            "lifecycle_event": self.lifecycle_event.as_dict(),
        }


@dataclass(frozen=True)
class ControlledBasisTrajectoryV253:
    initial_state: ThawedGaussianSpinorStateV252
    final_state: ThawedGaussianSpinorStateV252
    model: QuadraticSpinHamiltonianV252
    settings: ControlledBasisSettingsV253
    initial_packet_ids: tuple
    final_packet_ids: tuple
    initial_packet_ages: tuple
    final_packet_ages: tuple
    initial_next_packet_serial: int
    final_next_packet_serial: int
    steps: tuple
    claims: dict = field(default_factory=lambda: dict(V253_CONTROLLED_BASIS_CLAIMS))

    @property
    def event_counts(self):
        return {
            kind: sum(step.lifecycle_event.event_kind == kind for step in self.steps)
            for kind in ("none", "spawn", "prune", "merge")
        }

    @property
    def minimum_packet_count(self):
        return min(
            [self.initial_state.ngaussian] + [step.end.ngaussian for step in self.steps]
        )

    @property
    def maximum_packet_count(self):
        return max(
            [self.initial_state.ngaussian] + [step.end.ngaussian for step in self.steps]
        )

    @property
    def maximum_projection_loss(self):
        return float(max(
            (
                step.lifecycle_event.projection.relative_projection_loss
                for step in self.steps if step.lifecycle_event.projection is not None
            ),
            default=0.0,
        ))

    @property
    def maximum_event_energy_jump_hartree(self):
        return float(max(
            (
                abs(step.lifecycle_event.projection.energy_jump_hartree)
                for step in self.steps if step.lifecycle_event.projection is not None
            ),
            default=0.0,
        ))

    @property
    def maximum_nonlinear_residual(self):
        return float(max(
            (step.tdvp_step.nonlinear_residual_norm for step in self.steps),
            default=0.0,
        ))

    @property
    def maximum_norm_drift(self):
        initial = self.initial_state.generalized_norm
        values = [abs(step.tdvp_step.end.generalized_norm - initial) for step in self.steps]
        values.extend(abs(step.end.generalized_norm - initial) for step in self.steps)
        return float(max(values, default=0.0))

    def validate(self):
        self.model.validate()
        self.settings.validate()
        ids_initial, ages_initial = _validate_packet_metadata_v253(
            self.initial_packet_ids, self.initial_packet_ages, self.initial_state.ngaussian
        )
        ids_final, ages_final = _validate_packet_metadata_v253(
            self.final_packet_ids, self.final_packet_ages, self.final_state.ngaussian
        )
        serial_initial = _positive_integer_v253(
            self.initial_next_packet_serial, "initial_next_packet_serial", allow_zero=True
        )
        serial_final = _positive_integer_v253(
            self.final_next_packet_serial, "final_next_packet_serial", allow_zero=True
        )
        if type(self.claims) is not dict or any(type(value) is not bool for value in self.claims.values()):
            raise TypeError("every v0.25.3 claim must be a native Boolean.")
        if self.claims != V253_CONTROLLED_BASIS_CLAIMS:
            raise ValueError("v0.25.3 claims differ from the frozen boundary.")
        state = self.initial_state.validate()
        ids, ages, serial = ids_initial, ages_initial, serial_initial
        for expected_index, step in enumerate(self.steps, start=1):
            step.validate()
            if step.step_index != expected_index:
                raise ValueError("controlled trajectory step indices are not contiguous.")
            if _state_distance_v253(step.start, state) > self.settings.structural_tolerance:
                raise ValueError("controlled trajectory state chain is discontinuous.")
            if tuple(step.packet_ids_start) != ids or tuple(step.packet_ages_start) != ages or step.next_packet_serial_start != serial:
                raise ValueError("controlled trajectory metadata chain is discontinuous.")
            state = step.end
            ids, ages = tuple(step.packet_ids_end), tuple(step.packet_ages_end)
            serial = step.next_packet_serial_end
        if _state_distance_v253(state, self.final_state) > self.settings.structural_tolerance:
            raise ValueError("controlled trajectory final state is inconsistent.")
        if ids != ids_final or ages != ages_final or serial != serial_final:
            raise ValueError("controlled trajectory final metadata is inconsistent.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": CONTROLLED_BASIS_SCHEMA_V253,
            "residual_score": SPAWN_SCORE_V253,
            "projection_policy": PROJECTION_POLICY_V253,
            "event_order": EVENT_ORDER_V253,
            "potential_contract": POTENTIAL_CONTRACT_V253,
            "model": self.model.as_dict(),
            "settings": self.settings.as_dict(),
            "initial_state": self.initial_state.as_dict(),
            "final_state": self.final_state.as_dict(),
            "initial_packet_ids": list(self.initial_packet_ids),
            "final_packet_ids": list(self.final_packet_ids),
            "initial_packet_ages": list(self.initial_packet_ages),
            "final_packet_ages": list(self.final_packet_ages),
            "initial_next_packet_serial": int(self.initial_next_packet_serial),
            "final_next_packet_serial": int(self.final_next_packet_serial),
            "steps": [step.as_dict() for step in self.steps],
            "event_counts": self.event_counts,
            "minimum_packet_count": self.minimum_packet_count,
            "maximum_packet_count": self.maximum_packet_count,
            "maximum_projection_loss": self.maximum_projection_loss,
            "maximum_event_energy_jump_hartree": self.maximum_event_energy_jump_hartree,
            "maximum_nonlinear_residual": self.maximum_nonlinear_residual,
            "maximum_norm_drift": self.maximum_norm_drift,
            "claims": dict(self.claims),
        }

    def fingerprint(self):
        return _sha256_v252(self.as_dict())


def run_controlled_basis_dynamics_v253(
    initial_state,
    model,
    *,
    dt_au,
    steps,
    settings=ControlledBasisSettingsV253(),
):
    """Propagate v0.25.2 TDVP steps with controlled v0.25.3 topology events."""

    settings = settings.validate()
    model = model.validate()
    dt = float(dt_au)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("controlled basis forward dt must be finite and positive.")
    steps = _positive_integer_v253(steps, "controlled basis steps", allow_zero=True)
    state = initial_state.normalized()
    _validate_width_domain_v252(state, settings.tdvp_settings)
    if not settings.minimum_packet_count <= state.ngaussian <= settings.maximum_packet_count:
        raise ValueError("initial Gaussian count is outside the lifecycle domain.")
    if state.nstate != model.nstate:
        raise ValueError("controlled basis state/model electronic dimensions differ.")
    packet_ids = tuple(f"g{index:06d}" for index in range(state.ngaussian))
    packet_ages = tuple(0 for _ in range(state.ngaussian))
    next_packet_serial = state.ngaussian
    initial = state
    initial_ids = packet_ids
    initial_ages = packet_ages
    initial_serial = next_packet_serial
    receipts = []
    for step_index in range(1, steps + 1):
        start = state
        ids_start, ages_start, serial_start = packet_ids, packet_ages, next_packet_serial
        tdvp_step = controlled_tdvp_step_v253(
            state, model, dt, settings=settings
        )
        aged = tuple(value + 1 for value in packet_ages)
        if step_index % settings.adaptation_interval_steps == 0:
            event = adapt_basis_once_v253(
                tdvp_step.end,
                model,
                packet_ids=packet_ids,
                packet_ages=aged,
                next_packet_serial=next_packet_serial,
                settings=settings,
            )
        else:
            event = _no_event_v253(
                tdvp_step.end,
                model,
                settings,
                packet_ids,
                aged,
                next_packet_serial,
                reason="adaptation interval not reached",
            )
        state = event.after
        packet_ids = tuple(event.packet_ids_after)
        packet_ages = tuple(event.packet_ages_after)
        next_packet_serial = event.next_packet_serial_after
        receipts.append(
            ControlledBasisStepV253(
                step_index=step_index,
                start=start,
                end=state,
                packet_ids_start=ids_start,
                packet_ids_end=packet_ids,
                packet_ages_start=ages_start,
                packet_ages_end=packet_ages,
                next_packet_serial_start=serial_start,
                next_packet_serial_end=next_packet_serial,
                tdvp_step=tdvp_step,
                lifecycle_event=event,
                settings=settings,
            ).validate()
        )
    return ControlledBasisTrajectoryV253(
        initial_state=initial,
        final_state=state,
        model=model,
        settings=settings,
        initial_packet_ids=initial_ids,
        final_packet_ids=packet_ids,
        initial_packet_ages=initial_ages,
        final_packet_ages=packet_ages,
        initial_next_packet_serial=initial_serial,
        final_next_packet_serial=next_packet_serial,
        steps=tuple(receipts),
        claims=dict(V253_CONTROLLED_BASIS_CLAIMS),
    ).validate()
