"""Controlled multidimensional Gaussian-basis lifecycle for v0.26.0.

The event policy is intentionally conservative and deterministic:

    merge -> prune -> residual-driven spawn -> no event,

with at most one topology change at a checkpoint.  Spawn candidates are admitted
only after novelty, full-rank overlap, condition-number, and analytic TDVP-residual
gates.  Pruning and merging use fixed-time full-SVD projection with norm, fidelity,
projection-loss, and energy-jump receipts.

An exactly projected newborn has a zero coefficient and therefore null shape
tangents.  All coefficient coordinates remain active immediately, while its q, p,
log-width, and chirp coordinates remain frozen until its coefficient-row population
crosses the activation gate.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json

import numpy as np

from .multidimensional_gaussian_tdvp_v260 import (
    DiagonalGaussianSpinorStateV260,
    MultidimensionalImplicitMidpointStepV260,
    MultidimensionalVariationalSettingsV260,
    _cross_gaussian_data_v260,
    build_multidimensional_gaussian_matrices_v260,
    build_multidimensional_metric_system_v260,
    multidimensional_implicit_midpoint_step_v260,
    multidimensional_variational_energy_v260,
    pack_multidimensional_parameters_v260,
    residual_coupling_at_geometry_v260,
)


MULTIDIMENSIONAL_BASIS_SCHEMA_V260 = "gnd-controlled-multidimensional-basis-trajectory-v0.26.0"
MULTIDIMENSIONAL_EVENT_SCHEMA_V260 = "gnd-controlled-multidimensional-basis-event-v0.26.0"
SPAWN_SCORE_V260 = "norm of orthogonalized <g_candidate|dPsi/dt+iHPsi> divided by sqrt(novelty)"
PROJECTION_POLICY_V260 = "fixed-time full-SVD least-squares projection with fidelity, norm, loss, and energy receipts"
EVENT_ORDER_V260 = "merge then prune then spawn, with at most one event per checkpoint"
NEWBORN_ACTIVATION_V260 = (
    "coefficients active immediately; all shape coordinates dormant until both population "
    "and compatible-metric gates pass"
)


def _canonical_basis_v260(value):
    if isinstance(value, np.generic):
        return _canonical_basis_v260(value.item())
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.ndarray):
        return _canonical_basis_v260(value.tolist())
    if isinstance(value, dict):
        return {
            str(key): _canonical_basis_v260(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_basis_v260(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("basis-adaptation canonical data cannot be non-finite.")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported basis canonical value: {type(value).__name__}")


def _sha256_basis_v260(value):
    payload = json.dumps(
        _canonical_basis_v260(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _complex_pairs_basis_v260(value):
    array = np.asarray(value, dtype=complex)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def _positive_integer_basis_v260(value, name, *, allow_zero=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer.")
    minimum = 0 if allow_zero else 1
    if int(value) < minimum:
        raise ValueError(f"{name} must be at least {minimum}.")
    return int(value)


def _metadata_basis_v260(packet_ids, packet_ages, count):
    packet_ids = tuple(str(item) for item in packet_ids)
    packet_ages = tuple(_positive_integer_basis_v260(item, "packet age", allow_zero=True) for item in packet_ages)
    if len(packet_ids) != count or len(packet_ages) != count:
        raise ValueError("packet metadata length differs from the basis size.")
    if any(not item for item in packet_ids) or len(set(packet_ids)) != count:
        raise ValueError("packet IDs must be nonempty and unique.")
    return packet_ids, packet_ages


def _row_population_basis_v260(state, index):
    return float(np.sum(np.abs(state.coefficients[int(index)]) ** 2))


@dataclass(frozen=True)
class ControlledMultidimensionalBasisSettingsV260:
    tdvp_settings: MultidimensionalVariationalSettingsV260 = field(
        default_factory=MultidimensionalVariationalSettingsV260
    )
    basis_relative_cutoff: float = 1.0e-10
    basis_absolute_cutoff: float = 1.0e-12
    maximum_basis_condition_number: float = 1.0e10
    minimum_candidate_novelty: float = 1.0e-6
    spawn_residual_capture_threshold: float = 2.0e-5
    position_displacement_in_width_units: float = 1.0
    momentum_displacement_in_width_units: float = 1.0
    maximum_packet_count: int = 12
    # Inherited from the resolved v0.25.3 activation singularity: below this
    # population the shape-tangent block is numerically indistinguishable from
    # the metric null space and must remain dormant.
    newborn_activation_population: float = 1.0e-6
    maximum_activation_condition_number: float = 1.0e8
    maximum_activation_velocity_amplification: float = 100.0
    prune_population_threshold: float = 1.0e-10
    minimum_prune_age: int = 64
    minimum_merge_overlap: float = 0.997
    maximum_projection_loss: float = 2.0e-7
    maximum_event_energy_jump_hartree: float = 2.0e-6
    adapt_every_steps: int = 2
    spawn_enabled: bool = True
    prune_enabled: bool = True
    merge_enabled: bool = True
    maximum_events_per_checkpoint: int = 1
    event_order: str = EVENT_ORDER_V260
    spawn_score: str = SPAWN_SCORE_V260
    projection_policy: str = PROJECTION_POLICY_V260
    newborn_activation: str = NEWBORN_ACTIVATION_V260

    def validate(self):
        self.tdvp_settings.validate()
        positive = (
            "basis_relative_cutoff",
            "basis_absolute_cutoff",
            "maximum_basis_condition_number",
            "minimum_candidate_novelty",
            "spawn_residual_capture_threshold",
            "position_displacement_in_width_units",
            "momentum_displacement_in_width_units",
            "newborn_activation_population",
            "maximum_activation_condition_number",
            "maximum_activation_velocity_amplification",
            "prune_population_threshold",
            "minimum_merge_overlap",
            "maximum_projection_loss",
            "maximum_event_energy_jump_hartree",
        )
        for name in positive:
            value = float(getattr(self, name))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive.")
        if self.basis_relative_cutoff >= 1.0 or self.minimum_candidate_novelty >= 1.0:
            raise ValueError("relative cutoff and candidate novelty must be below one.")
        if self.minimum_merge_overlap > 1.0:
            raise ValueError("minimum merge overlap cannot exceed one.")
        _positive_integer_basis_v260(self.maximum_packet_count, "maximum_packet_count")
        _positive_integer_basis_v260(self.minimum_prune_age, "minimum_prune_age", allow_zero=True)
        _positive_integer_basis_v260(self.adapt_every_steps, "adapt_every_steps")
        _positive_integer_basis_v260(self.maximum_events_per_checkpoint, "maximum_events_per_checkpoint")
        for name in ("spawn_enabled", "prune_enabled", "merge_enabled"):
            if type(getattr(self, name)) is not bool:
                raise TypeError(f"{name} must be a native Boolean.")
            if not getattr(self, name):
                raise ValueError("v0.26.0 freezes spawn, prune, and merge as enabled.")
        if int(self.maximum_events_per_checkpoint) != 1:
            raise ValueError("v0.26.0 freezes at most one event per checkpoint.")
        if self.event_order != EVENT_ORDER_V260 or self.spawn_score != SPAWN_SCORE_V260:
            raise ValueError("v0.26.0 event order and spawn score are frozen.")
        if self.projection_policy != PROJECTION_POLICY_V260 or self.newborn_activation != NEWBORN_ACTIVATION_V260:
            raise ValueError("v0.26.0 projection and activation policies are frozen.")
        return self

    def as_dict(self):
        payload = asdict(self)
        payload["tdvp_settings"] = self.tdvp_settings.as_dict()
        return _canonical_basis_v260(payload)


@dataclass(frozen=True)
class MultidimensionalSpawnCandidateV260:
    q: np.ndarray
    p: np.ndarray
    widths: np.ndarray
    chirps: np.ndarray
    source_packet: int
    displacement_kind: str
    coordinate_axis: int
    sign: int

    def __post_init__(self):
        for name in ("q", "p", "widths", "chirps"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), dtype=float).copy())

    def validate(self, ndim=None):
        if self.q.ndim != 1 or len(self.q) < 1:
            raise ValueError("spawn candidate q must be a nonempty vector.")
        if any(value.shape != self.q.shape for value in (self.p, self.widths, self.chirps)):
            raise ValueError("spawn candidate geometry vectors differ in shape.")
        if ndim is not None and len(self.q) != int(ndim):
            raise ValueError("spawn candidate dimensionality differs from the state.")
        if not all(np.all(np.isfinite(value)) for value in (self.q, self.p, self.widths, self.chirps)):
            raise ValueError("spawn candidate contains non-finite values.")
        if np.min(self.widths) <= 0.0:
            raise ValueError("spawn candidate widths must be positive.")
        _positive_integer_basis_v260(self.source_packet, "source_packet", allow_zero=True)
        if self.displacement_kind not in ("position", "momentum"):
            raise ValueError("spawn displacement kind must be position or momentum.")
        if int(self.coordinate_axis) != self.coordinate_axis or not 0 <= int(self.coordinate_axis) < len(self.q):
            raise ValueError("spawn candidate coordinate axis is invalid.")
        if int(self.sign) not in (-1, 1):
            raise ValueError("spawn candidate sign must be -1 or +1.")
        return self

    def canonical_key(self):
        self.validate()
        return (
            int(self.source_packet),
            str(self.displacement_kind),
            int(self.coordinate_axis),
            int(self.sign),
            tuple(float(item) for item in self.q),
            tuple(float(item) for item in self.p),
            tuple(float(item) for item in self.widths),
            tuple(float(item) for item in self.chirps),
        )

    def as_dict(self):
        self.validate()
        return {
            "q": self.q.tolist(),
            "p": self.p.tolist(),
            "widths": self.widths.tolist(),
            "chirps": self.chirps.tolist(),
            "source_packet": int(self.source_packet),
            "displacement_kind": self.displacement_kind,
            "coordinate_axis": int(self.coordinate_axis),
            "sign": int(self.sign),
        }


def generate_multidimensional_spawn_candidates_v260(
    state, *, settings=ControlledMultidimensionalBasisSettingsV260()
):
    state = state.validate(require_normalized=False)
    settings = settings.validate()
    candidates = []
    for packet in range(state.ngaussian):
        for axis in range(state.ndim):
            position_step = float(settings.position_displacement_in_width_units) / np.sqrt(state.widths[packet, axis])
            momentum_step = float(settings.momentum_displacement_in_width_units) * np.sqrt(state.widths[packet, axis])
            for sign in (-1, 1):
                q = state.q[packet].copy()
                q[axis] += sign * position_step
                candidates.append(
                    MultidimensionalSpawnCandidateV260(
                        q, state.p[packet], state.widths[packet], state.chirps[packet],
                        packet, "position", axis, sign,
                    ).validate(state.ndim)
                )
                p = state.p[packet].copy()
                p[axis] += sign * momentum_step
                candidates.append(
                    MultidimensionalSpawnCandidateV260(
                        state.q[packet], p, state.widths[packet], state.chirps[packet],
                        packet, "momentum", axis, sign,
                    ).validate(state.ndim)
                )
    return tuple(sorted(candidates, key=lambda item: item.canonical_key()))


def _svd_basis_solve_v260(matrix, rhs, settings):
    matrix = np.asarray(matrix, dtype=complex)
    rhs = np.asarray(rhs, dtype=complex)
    left, singular_values, right_h = np.linalg.svd(matrix, full_matrices=True)
    maximum = float(singular_values[0])
    cutoff = max(float(settings.basis_absolute_cutoff), float(settings.basis_relative_cutoff) * maximum)
    retained = singular_values > cutoff
    rank = int(np.count_nonzero(retained))
    inverse = np.zeros_like(singular_values)
    inverse[retained] = 1.0 / singular_values[retained]
    solution = right_h.conj().T @ (inverse[:, None] * (left.conj().T @ np.atleast_2d(rhs).reshape(matrix.shape[0], -1)))
    if rhs.ndim == 1:
        solution = solution[:, 0]
    condition = float(maximum / singular_values[retained][-1]) if rank else float("inf")
    return solution, singular_values, cutoff, rank, condition


@dataclass(frozen=True)
class BasisProjectionReceiptV260:
    source_norm: float
    projected_norm_before_normalization: float
    source_projected_overlap: complex
    normalized_fidelity: float
    relative_projection_loss: float
    source_energy_hartree: float
    projected_energy_hartree: float
    target_overlap_singular_values: np.ndarray
    target_overlap_cutoff: float
    target_overlap_rank: int
    target_overlap_condition_number: float

    def __post_init__(self):
        object.__setattr__(self, "target_overlap_singular_values", np.asarray(self.target_overlap_singular_values, dtype=float).copy())

    @property
    def energy_jump_hartree(self):
        return float(self.projected_energy_hartree - self.source_energy_hartree)

    def validate(self):
        if self.target_overlap_singular_values.ndim != 1 or len(self.target_overlap_singular_values) < 1:
            raise ValueError("projection receipt requires a nonempty overlap spectrum.")
        if not np.all(np.isfinite(self.target_overlap_singular_values)) or np.min(self.target_overlap_singular_values) < 0.0:
            raise ValueError("projection overlap spectrum is invalid.")
        for name in (
            "source_norm", "projected_norm_before_normalization", "normalized_fidelity",
            "relative_projection_loss", "source_energy_hartree", "projected_energy_hartree",
            "target_overlap_cutoff", "target_overlap_condition_number",
        ):
            if not np.isfinite(float(getattr(self, name))):
                raise ValueError(f"projection receipt {name} is non-finite.")
        if self.source_norm <= 0.0 or self.projected_norm_before_normalization <= 0.0:
            raise ValueError("projection norms must be positive.")
        if not -2.0e-10 <= self.normalized_fidelity <= 1.0 + 2.0e-10:
            raise ValueError("projection fidelity is outside [0,1].")
        if self.relative_projection_loss < -2.0e-10:
            raise ValueError("projection loss cannot be negative.")
        expected_fidelity = abs(self.source_projected_overlap) ** 2 / (
            self.source_norm * self.projected_norm_before_normalization
        )
        expected_loss = (
            self.source_norm
            + self.projected_norm_before_normalization
            - 2.0 * float(np.real(self.source_projected_overlap))
        ) / self.source_norm
        if abs(self.normalized_fidelity - expected_fidelity) > 3.0e-10:
            raise ValueError("stored projection fidelity is inconsistent with its overlap receipt.")
        if abs(self.relative_projection_loss - max(expected_loss, 0.0)) > 3.0e-10:
            raise ValueError("stored relative projection loss is inconsistent with its overlap receipt.")
        if int(self.target_overlap_rank) != self.target_overlap_rank or not 1 <= int(self.target_overlap_rank) <= len(self.target_overlap_singular_values):
            raise ValueError("projection overlap rank is invalid.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "source_norm": float(self.source_norm),
            "projected_norm_before_normalization": float(self.projected_norm_before_normalization),
            "source_projected_overlap": [float(self.source_projected_overlap.real), float(self.source_projected_overlap.imag)],
            "normalized_fidelity": float(self.normalized_fidelity),
            "relative_projection_loss": float(self.relative_projection_loss),
            "source_energy_hartree": float(self.source_energy_hartree),
            "projected_energy_hartree": float(self.projected_energy_hartree),
            "energy_jump_hartree": self.energy_jump_hartree,
            "target_overlap_singular_values": self.target_overlap_singular_values.tolist(),
            "target_overlap_cutoff": float(self.target_overlap_cutoff),
            "target_overlap_rank": int(self.target_overlap_rank),
            "target_overlap_condition_number": float(self.target_overlap_condition_number),
        }


def project_multidimensional_state_v260(
    source,
    q,
    p,
    widths,
    chirps,
    model,
    *,
    settings=ControlledMultidimensionalBasisSettingsV260(),
):
    source = source.validate(require_normalized=False)
    model = model.validate()
    settings = settings.validate()
    q = np.asarray(q, dtype=float)
    p = np.asarray(p, dtype=float)
    widths = np.asarray(widths, dtype=float)
    chirps = np.asarray(chirps, dtype=float)
    if q.ndim != 2 or q.shape[0] < 1 or q.shape[1] != source.ndim:
        raise ValueError("target basis q has an invalid shape.")
    if any(value.shape != q.shape for value in (p, widths, chirps)) or np.min(widths) <= 0.0:
        raise ValueError("target basis packet arrays are incompatible.")
    target_count = q.shape[0]
    target_overlap = np.zeros((target_count, target_count), dtype=complex)
    cross = np.zeros((target_count, source.ngaussian), dtype=complex)
    for i in range(target_count):
        for j in range(target_count):
            target_overlap[i, j] = _cross_gaussian_data_v260(
                q[i], p[i], widths[i], chirps[i], q[j], p[j], widths[j], chirps[j]
            )[0]
        for j in range(source.ngaussian):
            cross[i, j] = _cross_gaussian_data_v260(
                q[i], p[i], widths[i], chirps[i],
                source.q[j], source.p[j], source.widths[j], source.chirps[j],
            )[0]
    rhs = cross @ source.coefficients
    coefficients, singular_values, cutoff, rank, condition = _svd_basis_solve_v260(target_overlap, rhs, settings)
    if rank != target_count:
        raise ValueError("target basis is rank deficient under the frozen SVD cutoff.")
    if condition > float(settings.maximum_basis_condition_number):
        raise ValueError("target basis exceeds the overlap condition-number gate.")
    projected_raw = DiagonalGaussianSpinorStateV260(
        q, p, widths, chirps, coefficients, source.time_au
    ).validate(require_normalized=False)
    source_norm = source.generalized_norm
    projected_norm = projected_raw.generalized_norm
    source_projected_overlap = complex(
        np.einsum("ia,ij,ja->", coefficients.conj(), cross, source.coefficients, optimize=True)
    )
    fidelity = float(abs(source_projected_overlap) ** 2 / (source_norm * projected_norm))
    residual_squared = max(
        source_norm + projected_norm - 2.0 * float(np.real(source_projected_overlap)), 0.0
    )
    loss = float(residual_squared / source_norm)
    projected = projected_raw.normalized()
    receipt = BasisProjectionReceiptV260(
        source_norm=source_norm,
        projected_norm_before_normalization=projected_norm,
        source_projected_overlap=source_projected_overlap,
        normalized_fidelity=min(max(fidelity, 0.0), 1.0),
        relative_projection_loss=loss,
        source_energy_hartree=multidimensional_variational_energy_v260(source, model),
        projected_energy_hartree=multidimensional_variational_energy_v260(projected, model),
        target_overlap_singular_values=singular_values,
        target_overlap_cutoff=cutoff,
        target_overlap_rank=rank,
        target_overlap_condition_number=condition,
    ).validate()
    return projected, receipt


@dataclass(frozen=True)
class SpawnCandidateEvaluationV260:
    candidate: MultidimensionalSpawnCandidateV260
    novelty: float
    residual_coupling: np.ndarray
    orthogonalized_residual_coupling: np.ndarray
    residual_capture: float
    enlarged_overlap_singular_values: np.ndarray
    enlarged_overlap_rank: int
    enlarged_overlap_condition_number: float
    rejection_reasons: tuple
    admitted: bool

    def __post_init__(self):
        object.__setattr__(self, "residual_coupling", np.asarray(self.residual_coupling, dtype=complex).copy())
        object.__setattr__(self, "orthogonalized_residual_coupling", np.asarray(self.orthogonalized_residual_coupling, dtype=complex).copy())
        object.__setattr__(self, "enlarged_overlap_singular_values", np.asarray(self.enlarged_overlap_singular_values, dtype=float).copy())
        object.__setattr__(self, "rejection_reasons", tuple(str(item) for item in self.rejection_reasons))

    def validate(self):
        self.candidate.validate()
        if not np.isfinite(float(self.novelty)) or self.novelty < -2.0e-10:
            raise ValueError("candidate novelty is invalid.")
        if self.residual_coupling.shape != self.orthogonalized_residual_coupling.shape or self.residual_coupling.ndim != 1:
            raise ValueError("candidate residual vectors have invalid shapes.")
        if not np.all(np.isfinite(self.residual_coupling)) or not np.all(np.isfinite(self.orthogonalized_residual_coupling)):
            raise ValueError("candidate residual vectors are non-finite.")
        if not np.isfinite(float(self.residual_capture)) or self.residual_capture < 0.0:
            raise ValueError("candidate residual capture is invalid.")
        if self.enlarged_overlap_singular_values.ndim != 1 or len(self.enlarged_overlap_singular_values) < 1:
            raise ValueError("candidate enlarged overlap spectrum is invalid.")
        if type(self.admitted) is not bool or self.admitted != (len(self.rejection_reasons) == 0):
            raise ValueError("candidate admission flag differs from its rejection reasons.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "candidate": self.candidate.as_dict(),
            "novelty": float(self.novelty),
            "residual_coupling": _complex_pairs_basis_v260(self.residual_coupling),
            "orthogonalized_residual_coupling": _complex_pairs_basis_v260(self.orthogonalized_residual_coupling),
            "residual_capture": float(self.residual_capture),
            "enlarged_overlap_singular_values": self.enlarged_overlap_singular_values.tolist(),
            "enlarged_overlap_rank": int(self.enlarged_overlap_rank),
            "enlarged_overlap_condition_number": float(self.enlarged_overlap_condition_number),
            "rejection_reasons": list(self.rejection_reasons),
            "admitted": self.admitted,
        }


def evaluate_multidimensional_spawn_candidate_v260(
    state,
    model,
    candidate,
    *,
    settings=ControlledMultidimensionalBasisSettingsV260(),
    metric_system=None,
):
    state = state.validate(require_normalized=False)
    model = model.validate()
    settings = settings.validate()
    candidate = candidate.validate(state.ndim)
    if candidate.source_packet >= state.ngaussian:
        raise ValueError("candidate source packet is outside the current basis.")
    if metric_system is None:
        active_mask = metric_compatible_activation_mask_v260(state, model, settings=settings)
        metric_system = build_multidimensional_metric_system_v260(
            state, model, settings=settings.tdvp_settings, active_shape_mask=active_mask
        )
    overlap = state.nuclear_overlap_matrix()
    b = np.asarray(
        [
            _cross_gaussian_data_v260(
                state.q[i], state.p[i], state.widths[i], state.chirps[i],
                candidate.q, candidate.p, candidate.widths, candidate.chirps,
            )[0]
            for i in range(state.ngaussian)
        ],
        dtype=complex,
    )
    coefficients, _, _, rank, condition = _svd_basis_solve_v260(overlap, b, settings)
    novelty = float(np.real(1.0 - np.vdot(b, coefficients)))
    novelty = max(novelty, 0.0)
    residual_candidate = residual_coupling_at_geometry_v260(
        state,
        model,
        metric_system.velocity,
        q=candidate.q,
        p=candidate.p,
        widths=candidate.widths,
        chirps=candidate.chirps,
    )
    basis_residual = np.asarray(
        [
            residual_coupling_at_geometry_v260(
                state,
                model,
                metric_system.velocity,
                q=state.q[i],
                p=state.p[i],
                widths=state.widths[i],
                chirps=state.chirps[i],
            )
            for i in range(state.ngaussian)
        ]
    )
    numerator = residual_candidate - np.einsum("i,ia->a", coefficients.conj(), basis_residual)
    if novelty > 0.0:
        orthogonalized = numerator / np.sqrt(novelty)
        residual_capture = float(np.linalg.norm(orthogonalized))
    else:
        orthogonalized = np.zeros_like(numerator)
        residual_capture = 0.0
    enlarged = np.block([[overlap, b[:, None]], [b.conj()[None, :], np.ones((1, 1), dtype=complex)]])
    singular_values = np.linalg.svd(enlarged, compute_uv=False)
    cutoff = max(settings.basis_absolute_cutoff, settings.basis_relative_cutoff * float(singular_values[0]))
    enlarged_rank = int(np.count_nonzero(singular_values > cutoff))
    enlarged_condition = (
        float(singular_values[0] / singular_values[enlarged_rank - 1]) if enlarged_rank else float("inf")
    )
    reasons = []
    if rank != state.ngaussian or condition > settings.maximum_basis_condition_number:
        reasons.append("current-basis-conditioning-failure")
    if novelty < settings.minimum_candidate_novelty:
        reasons.append("insufficient-novelty")
    if enlarged_rank != state.ngaussian + 1:
        reasons.append("rank-deficient-enlarged-basis")
    if enlarged_condition > settings.maximum_basis_condition_number:
        reasons.append("ill-conditioned-enlarged-basis")
    if residual_capture < settings.spawn_residual_capture_threshold:
        reasons.append("insufficient-residual-capture")
    return SpawnCandidateEvaluationV260(
        candidate=candidate,
        novelty=novelty,
        residual_coupling=residual_candidate,
        orthogonalized_residual_coupling=orthogonalized,
        residual_capture=residual_capture,
        enlarged_overlap_singular_values=singular_values,
        enlarged_overlap_rank=enlarged_rank,
        enlarged_overlap_condition_number=enlarged_condition,
        rejection_reasons=tuple(reasons),
        admitted=len(reasons) == 0,
    ).validate()


@dataclass(frozen=True)
class MultidimensionalBasisEventV260:
    event_kind: str
    reason: str
    before: DiagonalGaussianSpinorStateV260
    after: DiagonalGaussianSpinorStateV260
    packet_ids_before: tuple
    packet_ids_after: tuple
    packet_ages_before: tuple
    packet_ages_after: tuple
    next_packet_serial_before: int
    next_packet_serial_after: int
    projection: BasisProjectionReceiptV260 | None = None
    selected_candidate: MultidimensionalSpawnCandidateV260 | None = None
    candidate_evaluations: tuple = ()
    added_packet_id: str | None = None
    removed_packet_id: str | None = None

    def validate(self):
        if self.event_kind not in ("none", "spawn", "prune", "merge"):
            raise ValueError("basis event kind is invalid.")
        self.before.validate(require_normalized=False)
        self.after.validate(require_normalized=False)
        _metadata_basis_v260(self.packet_ids_before, self.packet_ages_before, self.before.ngaussian)
        _metadata_basis_v260(self.packet_ids_after, self.packet_ages_after, self.after.ngaussian)
        _positive_integer_basis_v260(self.next_packet_serial_before, "next serial", allow_zero=True)
        _positive_integer_basis_v260(self.next_packet_serial_after, "next serial", allow_zero=True)
        for item in self.candidate_evaluations:
            item.validate()
        if self.event_kind == "none":
            if self.before.ngaussian != self.after.ngaussian or self.projection is not None:
                raise ValueError("no-event receipt changed basis topology or stored a projection.")
            if np.max(
                np.abs(
                    pack_multidimensional_parameters_v260(self.before)
                    - pack_multidimensional_parameters_v260(self.after)
                )
            ) > 3.0e-10:
                raise ValueError("no-event receipt changed the variational state.")
        elif self.event_kind == "spawn":
            if self.after.ngaussian != self.before.ngaussian + 1 or self.selected_candidate is None:
                raise ValueError("spawn event did not add exactly one selected packet.")
            self.projection.validate()
            if self.added_packet_id is None or self.next_packet_serial_after != self.next_packet_serial_before + 1:
                raise ValueError("spawn event stable-ID metadata is incomplete.")
            if self.packet_ids_after[:-1] != self.packet_ids_before or self.packet_ids_after[-1] != self.added_packet_id:
                raise ValueError("spawn event packet-ID order is inconsistent.")
            if self.packet_ages_after[:-1] != self.packet_ages_before or self.packet_ages_after[-1] != 0:
                raise ValueError("spawn event packet ages are inconsistent.")
            if np.max(np.abs(self.after.coefficients[-1])) != 0.0:
                raise ValueError("an exactly projected newborn coefficient must be zero.")
            if not any(
                item.admitted
                and item.candidate.canonical_key() == self.selected_candidate.canonical_key()
                for item in self.candidate_evaluations
            ):
                raise ValueError("spawn selected candidate is not an admitted evaluation.")
        else:
            if self.after.ngaussian != self.before.ngaussian - 1 or self.removed_packet_id is None:
                raise ValueError("prune/merge event did not remove exactly one packet.")
            self.projection.validate()
            expected_ids = tuple(item for item in self.packet_ids_before if item != self.removed_packet_id)
            if self.packet_ids_after != expected_ids:
                raise ValueError("prune/merge packet IDs do not remove the declared packet.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIDIMENSIONAL_EVENT_SCHEMA_V260,
            "event_kind": self.event_kind,
            "reason": self.reason,
            "before": self.before.as_dict(),
            "after": self.after.as_dict(),
            "packet_ids_before": list(self.packet_ids_before),
            "packet_ids_after": list(self.packet_ids_after),
            "packet_ages_before": list(self.packet_ages_before),
            "packet_ages_after": list(self.packet_ages_after),
            "next_packet_serial_before": int(self.next_packet_serial_before),
            "next_packet_serial_after": int(self.next_packet_serial_after),
            "projection": None if self.projection is None else self.projection.as_dict(),
            "selected_candidate": None if self.selected_candidate is None else self.selected_candidate.as_dict(),
            "candidate_evaluations": [item.as_dict() for item in self.candidate_evaluations],
            "added_packet_id": self.added_packet_id,
            "removed_packet_id": self.removed_packet_id,
        }


def _identity_spawn_receipt_v260(before, after, model, settings):
    overlap = after.nuclear_overlap_matrix()
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    cutoff = max(settings.basis_absolute_cutoff, settings.basis_relative_cutoff * float(singular_values[0]))
    rank = int(np.count_nonzero(singular_values > cutoff))
    condition = float(singular_values[0] / singular_values[rank - 1])
    return BasisProjectionReceiptV260(
        source_norm=before.generalized_norm,
        projected_norm_before_normalization=after.generalized_norm,
        source_projected_overlap=complex(before.generalized_norm),
        normalized_fidelity=1.0,
        relative_projection_loss=0.0,
        source_energy_hartree=multidimensional_variational_energy_v260(before, model),
        projected_energy_hartree=multidimensional_variational_energy_v260(after, model),
        target_overlap_singular_values=singular_values,
        target_overlap_cutoff=cutoff,
        target_overlap_rank=rank,
        target_overlap_condition_number=condition,
    ).validate()


def metric_compatible_activation_mask_v260(
    state,
    model,
    *,
    settings=ControlledMultidimensionalBasisSettingsV260(),
    locked_active_mask=None,
):
    """Admit newborn shape blocks only when the SVD metric remains compatible.

    ``locked_active_mask`` preserves packets that were activated at an earlier
    step.  Dormant packets are tested one at a time in deterministic index order.
    A population threshold alone is insufficient in multiple dimensions because
    several shape directions can still lie below the retained metric scale.
    """

    state = state.validate(require_normalized=False)
    model = model.validate()
    settings = settings.validate()
    populations = np.sum(np.abs(state.coefficients) ** 2, axis=1)
    if locked_active_mask is None:
        mask = np.zeros(state.ngaussian, dtype=bool)
        mask[int(np.argmax(populations))] = True
    else:
        mask = np.asarray(locked_active_mask, dtype=bool).copy()
        if mask.shape != (state.ngaussian,):
            raise ValueError("locked active-shape mask has an invalid shape.")
        if not np.any(mask):
            mask[int(np.argmax(populations))] = True
    # Previously active directions are contractual; incompatibility here is a
    # real propagation failure and must not be hidden.
    base_system = build_multidimensional_metric_system_v260(
        state, model, settings=settings.tdvp_settings, active_shape_mask=mask
    )
    for packet in range(state.ngaussian):
        if mask[packet] or populations[packet] < settings.newborn_activation_population:
            continue
        trial = mask.copy()
        trial[packet] = True
        try:
            trial_system = build_multidimensional_metric_system_v260(
                state, model, settings=settings.tdvp_settings, active_shape_mask=trial
            )
        except ValueError:
            continue
        if (
            trial_system.solve_receipt.retained_condition_number
            > settings.maximum_activation_condition_number
        ):
            continue
        base_speed = max(base_system.solve_receipt.velocity_norm, 1.0e-14)
        if (
            trial_system.solve_receipt.velocity_norm
            > settings.maximum_activation_velocity_amplification * base_speed
        ):
            continue
        mask = trial
        base_system = trial_system
    return mask


def adapt_multidimensional_basis_once_v260(
    state,
    model,
    *,
    packet_ids=None,
    packet_ages=None,
    next_packet_serial=None,
    settings=ControlledMultidimensionalBasisSettingsV260(),
    active_shape_mask=None,
):
    state = state.validate(require_normalized=True, tolerance=settings.tdvp_settings.maximum_step_norm_drift)
    model = model.validate()
    settings = settings.validate()
    if packet_ids is None:
        packet_ids = tuple(f"g{index:06d}" for index in range(state.ngaussian))
    if packet_ages is None:
        packet_ages = tuple(0 for _ in range(state.ngaussian))
    packet_ids, packet_ages = _metadata_basis_v260(packet_ids, packet_ages, state.ngaussian)
    if next_packet_serial is None:
        next_packet_serial = state.ngaussian
    next_packet_serial = _positive_integer_basis_v260(next_packet_serial, "next_packet_serial", allow_zero=True)
    if f"g{next_packet_serial:06d}" in packet_ids:
        raise ValueError("next packet serial collides with an existing stable packet ID.")

    def no_event(reason, evaluations=()):
        return MultidimensionalBasisEventV260(
            "none", reason, state, state, packet_ids, packet_ids, packet_ages, packet_ages,
            next_packet_serial, next_packet_serial, candidate_evaluations=tuple(evaluations),
        ).validate()

    # Merge-to-survivor: deterministic pair order, remove the smaller coefficient row.
    overlap = state.nuclear_overlap_matrix()
    for i in range(state.ngaussian):
        for j in range(i + 1, state.ngaussian):
            if abs(overlap[i, j]) < settings.minimum_merge_overlap:
                continue
            remove = j if _row_population_basis_v260(state, j) <= _row_population_basis_v260(state, i) else i
            keep = np.asarray([index for index in range(state.ngaussian) if index != remove], dtype=int)
            projected, receipt = project_multidimensional_state_v260(
                state, state.q[keep], state.p[keep], state.widths[keep], state.chirps[keep],
                model, settings=settings,
            )
            if receipt.relative_projection_loss <= settings.maximum_projection_loss and abs(receipt.energy_jump_hartree) <= settings.maximum_event_energy_jump_hartree:
                return MultidimensionalBasisEventV260(
                    "merge", "overlap and projection gates passed", state, projected,
                    packet_ids, tuple(packet_ids[index] for index in keep),
                    packet_ages, tuple(packet_ages[index] for index in keep),
                    next_packet_serial, next_packet_serial,
                    projection=receipt, removed_packet_id=packet_ids[remove],
                ).validate()

    # Prune only an aged, coefficient-small packet that passes projection gates.
    prune_order = sorted(
        range(state.ngaussian), key=lambda index: (_row_population_basis_v260(state, index), packet_ids[index])
    )
    for remove in prune_order:
        if state.ngaussian <= 1:
            break
        if packet_ages[remove] < settings.minimum_prune_age or _row_population_basis_v260(state, remove) > settings.prune_population_threshold:
            continue
        keep = np.asarray([index for index in range(state.ngaussian) if index != remove], dtype=int)
        projected, receipt = project_multidimensional_state_v260(
            state, state.q[keep], state.p[keep], state.widths[keep], state.chirps[keep],
            model, settings=settings,
        )
        if receipt.relative_projection_loss <= settings.maximum_projection_loss and abs(receipt.energy_jump_hartree) <= settings.maximum_event_energy_jump_hartree:
            return MultidimensionalBasisEventV260(
                "prune", "population, age, and projection gates passed", state, projected,
                packet_ids, tuple(packet_ids[index] for index in keep),
                packet_ages, tuple(packet_ages[index] for index in keep),
                next_packet_serial, next_packet_serial,
                projection=receipt, removed_packet_id=packet_ids[remove],
            ).validate()

    if state.ngaussian >= settings.maximum_packet_count:
        return no_event("maximum packet count reached")
    if active_shape_mask is None:
        active_mask = metric_compatible_activation_mask_v260(state, model, settings=settings)
    else:
        active_mask = np.asarray(active_shape_mask, dtype=bool)
        if active_mask.shape != (state.ngaussian,):
            raise ValueError("active shape mask must contain one value per packet.")
    if np.any(~active_mask):
        # Serial activation is part of the controlled lifecycle: a second null
        # shape block is not introduced while an earlier newborn is still in
        # its coefficient-only stage.
        return no_event("dormant packet activation pending")
    metric_system = build_multidimensional_metric_system_v260(
        state, model, settings=settings.tdvp_settings, active_shape_mask=active_mask
    )
    evaluations = tuple(
        evaluate_multidimensional_spawn_candidate_v260(
            state, model, candidate, settings=settings, metric_system=metric_system
        )
        for candidate in generate_multidimensional_spawn_candidates_v260(state, settings=settings)
    )
    admitted = [item for item in evaluations if item.admitted]
    if not admitted:
        return no_event("no spawn candidate passed every gate", evaluations)
    selected = sorted(
        admitted, key=lambda item: (-item.residual_capture, item.candidate.canonical_key())
    )[0]
    candidate = selected.candidate
    after = DiagonalGaussianSpinorStateV260(
        np.vstack((state.q, candidate.q)),
        np.vstack((state.p, candidate.p)),
        np.vstack((state.widths, candidate.widths)),
        np.vstack((state.chirps, candidate.chirps)),
        np.vstack((state.coefficients, np.zeros((1, state.nstate), dtype=complex))),
        state.time_au,
    ).validate(require_normalized=True, tolerance=settings.tdvp_settings.maximum_step_norm_drift)
    added_id = f"g{next_packet_serial:06d}"
    receipt = _identity_spawn_receipt_v260(state, after, model, settings)
    return MultidimensionalBasisEventV260(
        "spawn", "highest residual candidate passed every gate", state, after,
        packet_ids, packet_ids + (added_id,), packet_ages, packet_ages + (0,),
        next_packet_serial, next_packet_serial + 1,
        projection=receipt, selected_candidate=candidate, candidate_evaluations=evaluations,
        added_packet_id=added_id,
    ).validate()


@dataclass(frozen=True)
class ControlledMultidimensionalStepV260:
    dynamics: MultidimensionalImplicitMidpointStepV260
    event: MultidimensionalBasisEventV260
    active_shape_mask: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, "active_shape_mask", np.asarray(self.active_shape_mask, dtype=bool).copy())

    def validate(self):
        self.dynamics.validate()
        self.event.validate()
        if self.active_shape_mask.shape != (self.dynamics.start.ngaussian,):
            raise ValueError("controlled-step active mask has an invalid shape.")
        if not np.array_equal(self.active_shape_mask, self.dynamics.active_shape_mask):
            raise ValueError("controlled-step active mask differs from the TDVP receipt.")
        if np.max(
            np.abs(
                pack_multidimensional_parameters_v260(self.dynamics.end)
                - pack_multidimensional_parameters_v260(self.event.before)
            )
        ) > self.dynamics.settings.structural_tolerance:
            raise ValueError("basis event does not start at the TDVP endpoint.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "dynamics": self.dynamics.as_dict(),
            "event": self.event.as_dict(),
            "active_shape_mask": self.active_shape_mask.tolist(),
        }


@dataclass(frozen=True)
class ControlledMultidimensionalTrajectoryV260:
    initial_state: DiagonalGaussianSpinorStateV260
    final_state: DiagonalGaussianSpinorStateV260
    model: object
    settings: ControlledMultidimensionalBasisSettingsV260
    steps: tuple
    initial_packet_ids: tuple
    final_packet_ids: tuple
    final_packet_ages: tuple
    next_packet_serial: int

    @property
    def event_counts(self):
        return {kind: sum(step.event.event_kind == kind for step in self.steps) for kind in ("none", "spawn", "prune", "merge")}

    @property
    def maximum_packet_count(self):
        return max([self.initial_state.ngaussian, *[step.event.after.ngaussian for step in self.steps]])

    @property
    def maximum_norm_drift(self):
        initial = self.initial_state.generalized_norm
        return float(max([abs(self.final_state.generalized_norm - initial), *[abs(step.event.after.generalized_norm - initial) for step in self.steps]], default=0.0))

    @property
    def maximum_projection_loss(self):
        return float(max([step.event.projection.relative_projection_loss for step in self.steps if step.event.projection is not None], default=0.0))

    def validate(self):
        settings = self.settings.validate()
        self.model.validate()
        self.initial_state.validate(require_normalized=True, tolerance=settings.tdvp_settings.maximum_step_norm_drift)
        self.final_state.validate(require_normalized=True, tolerance=settings.tdvp_settings.maximum_step_norm_drift)
        _metadata_basis_v260(self.initial_packet_ids, tuple(0 for _ in self.initial_packet_ids), self.initial_state.ngaussian)
        _metadata_basis_v260(self.final_packet_ids, self.final_packet_ages, self.final_state.ngaussian)
        current = self.initial_state
        packet_ids = self.initial_packet_ids
        packet_ages = tuple(0 for _ in packet_ids)
        for step in self.steps:
            step.validate()
            if current.q.shape != step.dynamics.start.q.shape or np.max(
                np.abs(
                    pack_multidimensional_parameters_v260(current)
                    - pack_multidimensional_parameters_v260(step.dynamics.start)
                )
            ) > settings.tdvp_settings.structural_tolerance:
                raise ValueError("controlled trajectory dynamics steps are not contiguous.")
            expected_ages = tuple(age + 1 for age in packet_ages)
            if step.event.packet_ids_before != packet_ids or step.event.packet_ages_before != expected_ages:
                raise ValueError("controlled trajectory packet metadata is not contiguous.")
            current = step.event.after
            packet_ids = step.event.packet_ids_after
            packet_ages = step.event.packet_ages_after
        if current.q.shape != self.final_state.q.shape or np.max(
            np.abs(
                pack_multidimensional_parameters_v260(current)
                - pack_multidimensional_parameters_v260(self.final_state)
            )
        ) > settings.tdvp_settings.structural_tolerance:
            raise ValueError("controlled trajectory final state is inconsistent.")
        if packet_ids != self.final_packet_ids or packet_ages != self.final_packet_ages:
            raise ValueError("controlled trajectory final packet metadata is inconsistent.")
        if self.maximum_packet_count > settings.maximum_packet_count:
            raise ValueError("controlled trajectory exceeded its packet cap.")
        return self

    def as_dict(self):
        self.validate()
        return {
            "schema": MULTIDIMENSIONAL_BASIS_SCHEMA_V260,
            "model_fingerprint": self.model.fingerprint(),
            "settings": self.settings.as_dict(),
            "initial_state": self.initial_state.as_dict(),
            "final_state": self.final_state.as_dict(),
            "steps": [step.as_dict() for step in self.steps],
            "initial_packet_ids": list(self.initial_packet_ids),
            "final_packet_ids": list(self.final_packet_ids),
            "final_packet_ages": list(self.final_packet_ages),
            "next_packet_serial": int(self.next_packet_serial),
            "event_counts": self.event_counts,
            "maximum_packet_count": self.maximum_packet_count,
            "maximum_norm_drift": self.maximum_norm_drift,
            "maximum_projection_loss": self.maximum_projection_loss,
        }

    def fingerprint(self):
        return _sha256_basis_v260(self.as_dict())


def run_controlled_multidimensional_dynamics_v260(
    initial_state,
    model,
    dt_au,
    steps,
    *,
    settings=ControlledMultidimensionalBasisSettingsV260(),
):
    settings = settings.validate()
    initial_state = initial_state.validate(require_normalized=True, tolerance=settings.tdvp_settings.maximum_step_norm_drift)
    steps = _positive_integer_basis_v260(steps, "steps", allow_zero=True)
    current = initial_state
    packet_ids = tuple(f"g{index:06d}" for index in range(current.ngaussian))
    packet_ages = tuple(0 for _ in range(current.ngaussian))
    initial_ids = packet_ids
    next_serial = current.ngaussian
    initial_populations = np.sum(np.abs(current.coefficients) ** 2, axis=1)
    active_ids = {packet_ids[int(np.argmax(initial_populations))]}
    receipts = []
    for step_index in range(1, steps + 1):
        locked_mask = np.asarray([packet_id in active_ids for packet_id in packet_ids], dtype=bool)
        active_mask = metric_compatible_activation_mask_v260(
            current,
            model,
            settings=settings,
            locked_active_mask=locked_mask,
        )
        active_ids.update(packet_id for packet_id, active in zip(packet_ids, active_mask) if active)
        dynamics = multidimensional_implicit_midpoint_step_v260(
            current,
            model,
            dt_au,
            settings=settings.tdvp_settings,
            active_shape_mask=active_mask,
        )
        aged = tuple(age + 1 for age in packet_ages)
        if step_index % settings.adapt_every_steps == 0:
            event = adapt_multidimensional_basis_once_v260(
                dynamics.end,
                model,
                packet_ids=packet_ids,
                packet_ages=aged,
                next_packet_serial=next_serial,
                settings=settings,
                active_shape_mask=active_mask,
            )
        else:
            event = MultidimensionalBasisEventV260(
                "none", "not an adaptation checkpoint", dynamics.end, dynamics.end,
                packet_ids, packet_ids, aged, aged, next_serial, next_serial,
            ).validate()
        receipts.append(ControlledMultidimensionalStepV260(dynamics, event, active_mask).validate())
        current = event.after
        packet_ids = event.packet_ids_after
        packet_ages = event.packet_ages_after
        next_serial = event.next_packet_serial_after
        active_ids.intersection_update(packet_ids)
    return ControlledMultidimensionalTrajectoryV260(
        initial_state=initial_state,
        final_state=current,
        model=model,
        settings=settings,
        steps=tuple(receipts),
        initial_packet_ids=initial_ids,
        final_packet_ids=packet_ids,
        final_packet_ages=packet_ages,
        next_packet_serial=next_serial,
    ).validate()


V260_MULTIDIMENSIONAL_BASIS_CLAIMS = {
    "multidimensional_residual_driven_spawning_validated": True,
    "multidimensional_projection_controlled_pruning_validated": True,
    "multidimensional_projection_controlled_merge_validated": True,
    "coefficient_only_newborn_activation_validated": True,
    "one_event_per_checkpoint_validated": True,
    "stable_packet_identity_and_age_validated": True,
    "packet_permutation_covariance_validated": True,
    "constant_electronic_gauge_covariance_validated": True,
    "arbitrary_multidimensional_spawning_directions_validated": False,
    "multiple_simultaneous_events_validated": False,
    "full_aims_branching_validated": False,
    "real_pyscf_soc_trajectory_admitted": False,
}
