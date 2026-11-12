"""Audited molecular electronic-frame bridge for the v0.30.0 development line.

This module binds the v0.29.3 coordinate-dependent geometry to an existing
``ElectronicOperatorSnapshotV21`` provider without weakening either contract.
The provider supplies total electronic Hamiltonians (including SOC when present),
physical operator derivatives, derivative connections/NACs, masses, and physical
cross-geometry overlaps.  Wilson transport is obtained either by integrating the
derivative connection or from the certified unitary polar factor of a physical
overlap.  The raw finite-manifold overlap is never used directly as a coefficient
transport.

The bridge is deliberately not a complete molecular trajectory propagator.  It is
the fail-closed interface needed before molecular data can enter the curved TDVP,
residual-spawning, and pre-ML data paths.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import json
from typing import Iterable

import numpy as np

from .curved_geometry_v293 import (
    ConnectionFieldV293,
    straight_path_v293,
    wilson_path_v293,
)
from .electronic_contract_v213 import ElectronicOperatorProvenanceV213
from .finite_manifold_transport_v233 import (
    CONSUMER_OVERLAP_POLICY_V233,
    FiniteManifoldOverlapPolicyV233,
    certify_reciprocal_transport_pair_v233,
)
from .pre_ml_contract_v294 import (
    SpawnOracleRecordV294,
    canonical_state_payload_v294,
    stable_hash_v294,
)
from .provider_numerical_identity_v233 import ProviderNumericalIdentityV233


V300_BRIDGE_SCHEMA = "gnd-molecular-frame-bridge-v0.30.0-dev"
V300_SEGMENT_SCHEMA = "gnd-molecular-frame-segment-audit-v1"
V300_CONTEXT_SCHEMA = "gnd-molecular-spawn-context-v1"


def _canonical_v300(value):
    if isinstance(value, np.generic):
        return _canonical_v300(value.item())
    if isinstance(value, np.ndarray):
        return _canonical_v300(value.tolist())
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("non-finite complex value in v0.30.0 payload")
        return [float(value.real), float(value.imag)]
    if isinstance(value, dict):
        return {
            str(key): _canonical_v300(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_v300(item) for item in value]
    if isinstance(value, float):
        if not np.isfinite(value):
            raise ValueError("non-finite real value in v0.30.0 payload")
        return float(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError(f"unsupported v0.30.0 payload value {type(value).__name__}")


def _hash_v300(value):
    raw = json.dumps(
        _canonical_v300(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _unitarity_residual_v300(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    identity = np.eye(matrix.shape[0], dtype=complex)
    return float(
        max(
            np.linalg.norm(matrix.conj().T @ matrix - identity, ord="fro"),
            np.linalg.norm(matrix @ matrix.conj().T - identity, ord="fro"),
        )
    )


def _resolve_provider_provenance_v300(provider):
    """Find inherited provenance through transparent provider wrappers."""
    seen = set()
    current = provider
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        provenance = getattr(current, "provenance", None)
        if provenance is not None:
            if type(provenance) is not ElectronicOperatorProvenanceV213:
                raise TypeError("provider provenance is not ElectronicOperatorProvenanceV213")
            return provenance.validate()
        current = getattr(current, "base_provider", None)
    raise ValueError("molecular frame provider lacks convention-complete provenance")


def _resolve_provider_identity_v300(provider, provenance, *, require_numerical_identity):
    # Numerical identity is not inherited through wrappers.  A gauge, tracking, or
    # interpolation wrapper changes numerical behavior and must explicitly attest its
    # own complete identity even when physical provenance is inherited from its base.
    identity = getattr(provider, "numerical_identity_v233", None)
    identity = identity() if callable(identity) else identity
    if identity is not None:
        if type(identity) is not ProviderNumericalIdentityV233:
            raise TypeError("provider numerical identity has the wrong type")
        identity = identity.validate()
        if identity.provenance_fingerprint != provenance.fingerprint():
            raise ValueError("provider numerical identity and provenance disagree")
        return "v0.23.3-numerical-identity", identity.fingerprint()
    if require_numerical_identity:
        raise ValueError(
            "provider lacks v0.23.3 numerical identity; production molecular "
            "bridges fail closed"
        )
    return "provenance-only-development", provenance.fingerprint()


@dataclass(frozen=True)
class MolecularFrameBridgeSettingsV300:
    coordinate_tolerance: float = 1.0e-11
    structural_tolerance: float = 1.0e-12
    reciprocity_tolerance: float = 2.0e-9
    connection_overlap_tolerance: float = 2.0e-7
    maximum_audit_segment_length: float = 0.5
    require_numerical_identity: bool = False
    overlap_policy: FiniteManifoldOverlapPolicyV233 = field(
        default_factory=lambda: CONSUMER_OVERLAP_POLICY_V233
    )

    def validate(self):
        for name in (
            "coordinate_tolerance",
            "structural_tolerance",
            "reciprocity_tolerance",
            "connection_overlap_tolerance",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        length = float(self.maximum_audit_segment_length)
        if not np.isfinite(length) or length <= 0.0:
            raise ValueError("maximum_audit_segment_length must be finite and positive")
        if not isinstance(self.require_numerical_identity, (bool, np.bool_)):
            raise TypeError("require_numerical_identity must be Boolean")
        self.overlap_policy.validate()
        return self


@dataclass(frozen=True)
class MolecularFrameSegmentAuditV300:
    schema: str
    start: np.ndarray
    end: np.ndarray
    segment_length: float
    overlap_singular_values: np.ndarray
    overlap_condition_number: float
    maximum_principal_angle_radians: float
    overlap_reciprocity_residual: float
    transport_reciprocity_residual: float
    overlap_transport_unitarity_residual: float
    connection_transport_unitarity_residual: float
    connection_overlap_residual: float
    provider_fingerprint: str
    bridge_fingerprint: str
    passed: bool

    def validate(self):
        if self.schema != V300_SEGMENT_SCHEMA:
            raise ValueError("molecular segment audit schema mismatch")
        start = np.asarray(self.start, dtype=float)
        end = np.asarray(self.end, dtype=float)
        singular = np.asarray(self.overlap_singular_values, dtype=float)
        if start.ndim != 1 or end.shape != start.shape:
            raise ValueError("segment endpoints have incompatible shapes")
        if singular.ndim != 1 or singular.size < 1:
            raise ValueError("segment audit requires overlap singular values")
        values = (
            self.segment_length,
            self.overlap_condition_number,
            self.maximum_principal_angle_radians,
            self.overlap_reciprocity_residual,
            self.transport_reciprocity_residual,
            self.overlap_transport_unitarity_residual,
            self.connection_transport_unitarity_residual,
            self.connection_overlap_residual,
        )
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(singular)):
            raise ValueError("segment audit contains non-finite diagnostics")
        for fingerprint in (self.provider_fingerprint, self.bridge_fingerprint):
            if len(fingerprint) != 64:
                raise ValueError("segment audit fingerprints must be SHA-256 digests")
        if not isinstance(self.passed, (bool, np.bool_)):
            raise TypeError("segment audit decision must be Boolean")
        return self

    def as_dict(self):
        self.validate()
        return _canonical_v300(asdict(self))


@dataclass(frozen=True)
class MolecularSpawnContextV300:
    schema: str
    base_oracle_record_fingerprint: str
    provider_fingerprint: str
    bridge_fingerprint: str
    provider_identity_kind: str
    path_family: str
    packet_ids: tuple
    point_fingerprints: tuple
    adiabatic_energy_spectra: tuple
    minimum_adiabatic_gaps: tuple
    total_operator_gradient_norms: tuple
    soc_enabled: bool
    state_fingerprint: str
    trajectory_group_id: str
    split: str

    def validate(self):
        if self.schema != V300_CONTEXT_SCHEMA:
            raise ValueError("molecular spawn context schema mismatch")
        n = len(self.packet_ids)
        if n < 1 or len(set(self.packet_ids)) != n:
            raise ValueError("molecular context packet IDs must be nonempty and unique")
        for values in (
            self.point_fingerprints,
            self.adiabatic_energy_spectra,
            self.minimum_adiabatic_gaps,
            self.total_operator_gradient_norms,
        ):
            if len(values) != n:
                raise ValueError("molecular context packet arrays have different lengths")
        for fingerprint in (
            self.base_oracle_record_fingerprint,
            self.provider_fingerprint,
            self.bridge_fingerprint,
            self.state_fingerprint,
            *self.point_fingerprints,
        ):
            if len(str(fingerprint)) != 64:
                raise ValueError("molecular context fingerprints must be SHA-256 digests")
        if self.split not in ("train", "validation", "test"):
            raise ValueError("molecular context has an invalid dataset split")
        if not isinstance(self.soc_enabled, (bool, np.bool_)):
            raise TypeError("soc_enabled must be Boolean")
        _canonical_v300(self.as_dict_unchecked())
        return self

    def as_dict_unchecked(self):
        return asdict(self)

    def as_dict(self):
        return _canonical_v300(self.as_dict_unchecked())

    @property
    def fingerprint(self):
        self.validate()
        return _hash_v300(self.as_dict())


class MolecularFrameBridgeV300:
    """Validate and expose one fixed-dimension molecular electronic provider."""

    def __init__(
        self,
        provider,
        reference_q,
        *,
        settings=MolecularFrameBridgeSettingsV300(),
    ):
        self.provider = provider
        self.settings = settings.validate()
        self.provenance = _resolve_provider_provenance_v300(provider)
        self.provider_identity_kind, self.provider_fingerprint = (
            _resolve_provider_identity_v300(
                provider,
                self.provenance,
                require_numerical_identity=self.settings.require_numerical_identity,
            )
        )
        reference_q = np.asarray(reference_q, dtype=float)
        if reference_q.ndim != 1 or reference_q.size < 1 or not np.all(np.isfinite(reference_q)):
            raise ValueError("reference_q must be a finite nonempty coordinate vector")
        self.reference_q = reference_q.copy()
        reference = self._raw_snapshot(reference_q)
        self.ndim = reference.point.nq
        self.nstate = reference.point.nstate
        if self.provenance.model_space.nstate != self.nstate:
            raise ValueError("provider point and provenance model-space dimensions differ")
        self._reference_mass = reference.point.mass_matrix_q_au.copy()
        self.bridge_fingerprint = _hash_v300(
            {
                "schema": V300_BRIDGE_SCHEMA,
                "provider_class": (
                    f"{provider.__class__.__module__}.{provider.__class__.__qualname__}"
                ),
                "provider_identity_kind": self.provider_identity_kind,
                "provider_fingerprint": self.provider_fingerprint,
                "reference_q": self.reference_q,
                "reference_point": self._point_payload(reference.point),
                "settings": {
                    **asdict(self.settings),
                    "overlap_policy": self.settings.overlap_policy.as_dict(),
                },
            }
        )

    def _raw_snapshot(self, q):
        if not hasattr(self.provider, "evaluate_snapshot"):
            raise TypeError("molecular frame provider must implement evaluate_snapshot(q)")
        snapshot = self.provider.evaluate_snapshot(np.asarray(q, dtype=float))
        if not hasattr(snapshot, "validate") or not hasattr(snapshot, "point"):
            raise TypeError("provider did not return an electronic operator snapshot")
        return snapshot.validate(atol=self.settings.structural_tolerance)

    @staticmethod
    def _point_payload(point):
        return {
            "q": point.q,
            "H": point.H,
            "dH_operator_dq": point.hamiltonian_derivative_operator_q,
            "connection_q": point.connection_q,
            "mass_matrix_q_au": point.mass_matrix_q_au,
        }

    def evaluate_snapshot(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (self.ndim,) or not np.all(np.isfinite(q)):
            raise ValueError("molecular bridge coordinate has the wrong shape or is non-finite")
        snapshot = self._raw_snapshot(q)
        point = snapshot.point
        if point.nq != self.ndim or point.nstate != self.nstate:
            raise ValueError("molecular provider changed coordinate or electronic dimension")
        if np.linalg.norm(point.q - q) > self.settings.coordinate_tolerance:
            raise ValueError("molecular provider returned a point at the wrong coordinate")
        if point.mass_matrix_q_au.shape != self._reference_mass.shape:
            raise ValueError("molecular provider changed generalized mass dimension")
        return snapshot

    def operator_point(self, q):
        return self.evaluate_snapshot(q).point

    def parameter_jacobian(self, state, connection, quadrature, *, path_factory, active_shape_mask):
        """Optional weighted parameter Jacobian; None uses central differences.

        A specialized bridge must validate this derivative independently and bind
        its backend convention into its bridge fingerprint. Generic molecular
        providers retain the existing finite-difference tangent.
        """
        return None

    def point_fingerprint(self, q):
        point = self.operator_point(q)
        return _hash_v300(
            {
                "schema": V300_BRIDGE_SCHEMA,
                "provider_fingerprint": self.provider_fingerprint,
                "point": self._point_payload(point),
            }
        )

    def connection_field(self, *, transport_mode="nac_ode"):
        """Return the v0.29.3 connection field in one explicit transport mode.

        ``nac_ode`` integrates the provider derivative connection along the path.
        ``overlap_polar`` uses the certified polar factor for every path segment.
        """
        if transport_mode not in {"nac_ode", "overlap_polar"}:
            raise ValueError("transport_mode must be 'nac_ode' or 'overlap_polar'")

        def evaluator(q):
            return self.operator_point(q).connection_q

        segment_transport = None
        if transport_mode == "overlap_polar":
            segment_transport = lambda start, end: self.overlap_transport(start, end)
        return ConnectionFieldV293(
            self.ndim,
            self.nstate,
            evaluator,
            name=f"molecular-v300-{transport_mode}-{self.bridge_fingerprint[:12]}",
            segment_transport=segment_transport,
        ).validate([self.reference_q])

    def _reciprocal_transport(self, start, end):
        start_snapshot = self.evaluate_snapshot(start)
        end_snapshot = self.evaluate_snapshot(end)
        if not hasattr(self.provider, "snapshot_overlap"):
            raise TypeError("molecular frame provider must implement snapshot_overlap")
        # <Phi(end)|Phi(start)> maps coefficients at start into the end frame after
        # taking its certified unitary polar factor.
        end_start = np.asarray(
            self.provider.snapshot_overlap(end_snapshot, start_snapshot), dtype=complex
        )
        start_end = np.asarray(
            self.provider.snapshot_overlap(start_snapshot, end_snapshot), dtype=complex
        )
        return certify_reciprocal_transport_pair_v233(
            end_start,
            start_end,
            policy=self.settings.overlap_policy,
            reciprocity_tolerance=self.settings.reciprocity_tolerance,
        )

    def overlap_transport(self, start, end):
        pair = self._reciprocal_transport(start, end)
        return pair.left_to_right_block.right_to_left_transport.copy()

    def audit_segment(self, start, end, *, require_consistency=True):
        start = np.asarray(start, dtype=float)
        end = np.asarray(end, dtype=float)
        if start.shape != (self.ndim,) or end.shape != (self.ndim,):
            raise ValueError("molecular audit endpoints have the wrong shape")
        length = float(np.linalg.norm(end - start))
        if length <= 0.0:
            raise ValueError("molecular audit segment must have nonzero length")
        if length > self.settings.maximum_audit_segment_length:
            raise ValueError("molecular audit segment exceeds the configured local limit")
        pair = self._reciprocal_transport(start, end)
        overlap_result = pair.left_to_right_block
        overlap_transport = overlap_result.right_to_left_transport
        connection_transport = wilson_path_v293(
            self.connection_field(transport_mode="nac_ode"),
            straight_path_v293(start, end),
        )
        residual = float(
            np.linalg.norm(connection_transport - overlap_transport, ord="fro")
            / max(np.linalg.norm(overlap_transport, ord="fro"), 1.0)
        )
        passed = bool(residual <= self.settings.connection_overlap_tolerance)
        receipt = MolecularFrameSegmentAuditV300(
            schema=V300_SEGMENT_SCHEMA,
            start=start.copy(),
            end=end.copy(),
            segment_length=length,
            overlap_singular_values=overlap_result.singular_values.copy(),
            overlap_condition_number=float(overlap_result.condition_number),
            maximum_principal_angle_radians=float(
                overlap_result.maximum_principal_angle_radians
            ),
            overlap_reciprocity_residual=float(pair.overlap_reciprocity_residual),
            transport_reciprocity_residual=float(pair.transport_reciprocity_residual),
            overlap_transport_unitarity_residual=_unitarity_residual_v300(
                overlap_transport
            ),
            connection_transport_unitarity_residual=_unitarity_residual_v300(
                connection_transport
            ),
            connection_overlap_residual=residual,
            provider_fingerprint=self.provider_fingerprint,
            bridge_fingerprint=self.bridge_fingerprint,
            passed=passed,
        ).validate()
        if require_consistency and not passed:
            raise ValueError(
                "provider derivative connection and overlap transport disagree on "
                "the audited segment"
            )
        return receipt

    def audit_path(self, vertices: Iterable[np.ndarray], *, require_consistency=True):
        vertices = [np.asarray(value, dtype=float) for value in vertices]
        if len(vertices) < 2:
            raise ValueError("molecular path audit requires at least two vertices")
        return tuple(
            self.audit_segment(left, right, require_consistency=require_consistency)
            for left, right in zip(vertices[:-1], vertices[1:])
        )

    def apply_electronic_hamiltonian(self, points, psi):
        """Apply only the molecular electronic Hamiltonian pointwise.

        The nuclear kinetic contribution remains the responsibility of the inherited
        grid/TDVP Hamiltonian.  This explicit name prevents an electronic-only action
        from being mistaken for the full molecular Hamiltonian.
        """
        points = np.asarray(points, dtype=float)
        psi = np.asarray(psi, dtype=complex)
        if points.ndim != 2 or points.shape[1] != self.ndim:
            raise ValueError("electronic-action points have the wrong shape")
        if psi.shape != (points.shape[0], self.nstate):
            raise ValueError("electronic-action wavefunction has the wrong shape")
        out = np.empty_like(psi)
        for index, point in enumerate(points):
            out[index] = self.operator_point(point).H @ psi[index]
        return out


def molecular_spawn_context_v300(
    record: SpawnOracleRecordV294,
    state,
    bridge: MolecularFrameBridgeV300,
    *,
    packet_ids=None,
    path_family="straight_path_v293",
):
    """Bind a frozen v0.29.4 oracle record to molecular provider provenance.

    Only gauge-invariant spectra and operator-gradient norms are exposed as model
    context.  Raw Hamiltonian, NAC, SOC, and coefficient components remain excluded.
    """
    if type(record) is not SpawnOracleRecordV294:
        raise TypeError("record must be SpawnOracleRecordV294")
    record = record.validate()
    state = state.validate(require_normalized=False)
    if (state.ndim, state.nstate) != (bridge.ndim, bridge.nstate):
        raise ValueError("state and molecular bridge dimensions disagree")
    ids = tuple(packet_ids or (f"g{index:06d}" for index in range(state.ngaussian)))
    expected_state = stable_hash_v294(
        canonical_state_payload_v294(
            state, packet_ids=ids, path_family=str(path_family)
        )
    )
    if record.state_fingerprint != expected_state:
        raise ValueError("molecular context state does not match the oracle record")
    order = np.argsort(np.asarray(ids, dtype=str), kind="stable")
    ordered_ids = tuple(ids[index] for index in order)
    point_fingerprints = []
    spectra = []
    gaps = []
    gradient_norms = []
    for index in order:
        point = bridge.operator_point(state.q[index])
        energies = np.sort(np.linalg.eigvalsh(point.H))
        point_fingerprints.append(bridge.point_fingerprint(state.q[index]))
        spectra.append(tuple(float(value) for value in energies))
        gaps.append(
            float(np.min(np.diff(energies))) if len(energies) > 1 else 0.0
        )
        gradient_norms.append(float(np.linalg.norm(point.dH_dq)))
    context = MolecularSpawnContextV300(
        schema=V300_CONTEXT_SCHEMA,
        base_oracle_record_fingerprint=record.record_fingerprint,
        provider_fingerprint=bridge.provider_fingerprint,
        bridge_fingerprint=bridge.bridge_fingerprint,
        provider_identity_kind=bridge.provider_identity_kind,
        path_family=str(path_family),
        packet_ids=ordered_ids,
        point_fingerprints=tuple(point_fingerprints),
        adiabatic_energy_spectra=tuple(spectra),
        minimum_adiabatic_gaps=tuple(gaps),
        total_operator_gradient_norms=tuple(gradient_norms),
        soc_enabled=bool(bridge.provenance.soc_enabled),
        state_fingerprint=record.state_fingerprint,
        trajectory_group_id=record.trajectory_group_id,
        split=record.split,
    )
    return context.validate()


def serialize_molecular_spawn_contexts_v300(contexts):
    contexts = [context.validate() for context in contexts]
    rows = sorted(
        (context.as_dict() for context in contexts),
        key=lambda row: (
            row["trajectory_group_id"],
            row["state_fingerprint"],
            row["base_oracle_record_fingerprint"],
        ),
    )
    return "".join(
        json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n"
        for row in rows
    )
