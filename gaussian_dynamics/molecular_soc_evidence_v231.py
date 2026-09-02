"""Raw-observation molecular SOC evidence for v0.23.1.

v0.23.0 froze the evidence categories.  v0.23.1 derives their numerical values
from stored observations instead of accepting precomputed pass/fail assertions.
"""

from dataclasses import dataclass
import numpy as np

from .molecular_soc_contract_v230 import MolecularSOCValidationEvidenceV230


_UNITS_V231 = {"hartree", "hartree/bohr", "dimensionless"}
_METRICS_V231 = {"maximum_absolute", "root_mean_square"}


def _required_text_v231(name, value):
    text = str(value).strip()
    if not text:
        raise ValueError(f"{name} cannot be empty.")
    return text


def _finite_nonnegative_v231(name, value):
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative.")
    return result


def _probability_v231(name, value):
    result = float(value)
    if not np.isfinite(result) or result < 0.0 or result > 1.0:
        raise ValueError(f"{name} must lie in [0,1].")
    return result


def _native_bool_v231(name, value):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be Boolean.")
    return bool(value)


def _complex_tuple_v231(name, values, shape):
    shape = tuple(int(value) for value in shape)
    if not shape or any(value < 1 for value in shape):
        raise ValueError(f"{name}_shape must contain positive dimensions.")
    array = np.asarray(values, dtype=complex)
    if array.size != int(np.prod(shape)):
        raise ValueError(f"{name} size disagrees with its declared shape.")
    array = array.reshape(shape)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values.")
    return array


def _metric_v231(left, right, metric):
    difference = np.asarray(left, dtype=complex) - np.asarray(right, dtype=complex)
    if metric == "maximum_absolute":
        return float(np.max(np.abs(difference), initial=0.0))
    if metric == "root_mean_square":
        return float(np.sqrt(np.mean(np.abs(difference) ** 2)))
    raise ValueError(f"unsupported evidence metric {metric!r}.")


def _complex_list_v231(values):
    return [
        {"real": float(complex(value).real), "imag": float(complex(value).imag)}
        for value in np.asarray(values, dtype=complex).reshape(-1)
    ]


def _complex_from_list_v231(values):
    return tuple(complex(float(value["real"]), float(value["imag"])) for value in values)


@dataclass(frozen=True)
class IndependentReferenceObservationV231:
    reference_id: str
    observable: str
    unit: str
    value_shape: tuple[int, ...]
    computed_values: tuple[complex, ...]
    reference_values: tuple[complex, ...]
    computed_artifact: str
    reference_artifact: str
    tolerance: float
    metric: str = "maximum_absolute"

    def validate(self):
        _required_text_v231("reference_id", self.reference_id)
        _required_text_v231("observable", self.observable)
        _required_text_v231("computed_artifact", self.computed_artifact)
        _required_text_v231("reference_artifact", self.reference_artifact)
        if self.computed_artifact == self.reference_artifact:
            raise ValueError("computed and independent-reference artifacts must differ.")
        if self.unit not in _UNITS_V231:
            raise ValueError(f"unsupported evidence unit {self.unit!r}.")
        if self.metric not in _METRICS_V231:
            raise ValueError(f"unsupported evidence metric {self.metric!r}.")
        _finite_nonnegative_v231("reference tolerance", self.tolerance)
        _complex_tuple_v231("computed_values", self.computed_values, self.value_shape)
        _complex_tuple_v231("reference_values", self.reference_values, self.value_shape)
        return self

    @property
    def error(self):
        self.validate()
        computed = _complex_tuple_v231(
            "computed_values", self.computed_values, self.value_shape
        )
        reference = _complex_tuple_v231(
            "reference_values", self.reference_values, self.value_shape
        )
        return _metric_v231(computed, reference, self.metric)

    @property
    def passed(self):
        return bool(self.error <= float(self.tolerance))

    def as_dict(self):
        self.validate()
        return {
            "reference_id": self.reference_id,
            "observable": self.observable,
            "unit": self.unit,
            "value_shape": list(self.value_shape),
            "computed_values": _complex_list_v231(self.computed_values),
            "reference_values": _complex_list_v231(self.reference_values),
            "computed_artifact": self.computed_artifact,
            "reference_artifact": self.reference_artifact,
            "tolerance": float(self.tolerance),
            "metric": self.metric,
            "derived_error": self.error,
            "passed": self.passed,
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            reference_id=payload["reference_id"],
            observable=payload["observable"],
            unit=payload["unit"],
            value_shape=tuple(payload["value_shape"]),
            computed_values=_complex_from_list_v231(payload["computed_values"]),
            reference_values=_complex_from_list_v231(payload["reference_values"]),
            computed_artifact=payload["computed_artifact"],
            reference_artifact=payload["reference_artifact"],
            tolerance=payload["tolerance"],
            metric=payload["metric"],
        ).validate()


@dataclass(frozen=True)
class ConvergenceLadderObservationV231:
    kind: str
    labels: tuple[str, ...]
    observable: str
    unit: str
    value_shape: tuple[int, ...]
    values: tuple[tuple[complex, ...], ...]
    source_artifacts: tuple[str, ...]
    tolerance: float
    metric: str = "maximum_absolute"

    def validate(self):
        if self.kind not in {"basis", "method"}:
            raise ValueError("convergence kind must be 'basis' or 'method'.")
        labels = tuple(str(value).strip() for value in self.labels)
        if len(labels) < 2 or any(not value for value in labels):
            raise ValueError("a convergence ladder requires at least two labels.")
        if len(set(labels)) != len(labels):
            raise ValueError("convergence ladder labels must be unique.")
        if len(self.values) != len(labels):
            raise ValueError("convergence ladder requires one value per label.")
        if len(self.source_artifacts) != len(labels):
            raise ValueError("convergence ladder requires one source artifact per label.")
        if len(set(self.source_artifacts)) != len(self.source_artifacts):
            raise ValueError("convergence source artifacts must be distinct.")
        for artifact in self.source_artifacts:
            _required_text_v231("convergence source artifact", artifact)
        _required_text_v231("convergence observable", self.observable)
        if self.unit not in _UNITS_V231:
            raise ValueError(f"unsupported evidence unit {self.unit!r}.")
        if self.metric not in _METRICS_V231:
            raise ValueError(f"unsupported evidence metric {self.metric!r}.")
        _finite_nonnegative_v231("convergence tolerance", self.tolerance)
        for index, value in enumerate(self.values):
            _complex_tuple_v231(f"convergence value {index}", value, self.value_shape)
        return self

    @property
    def changes(self):
        self.validate()
        arrays = [
            _complex_tuple_v231("convergence value", value, self.value_shape)
            for value in self.values
        ]
        return tuple(
            _metric_v231(left, right, self.metric)
            for left, right in zip(arrays[:-1], arrays[1:])
        )

    @property
    def passed(self):
        return bool(self.changes[-1] <= float(self.tolerance))

    def as_dict(self):
        self.validate()
        return {
            "kind": self.kind,
            "labels": list(self.labels),
            "observable": self.observable,
            "unit": self.unit,
            "value_shape": list(self.value_shape),
            "values": [_complex_list_v231(value) for value in self.values],
            "source_artifacts": list(self.source_artifacts),
            "tolerance": float(self.tolerance),
            "metric": self.metric,
            "derived_changes": list(self.changes),
            "passed": self.passed,
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            kind=payload["kind"],
            labels=tuple(payload["labels"]),
            observable=payload["observable"],
            unit=payload["unit"],
            value_shape=tuple(payload["value_shape"]),
            values=tuple(
                _complex_from_list_v231(value) for value in payload["values"]
            ),
            source_artifacts=tuple(payload["source_artifacts"]),
            tolerance=payload["tolerance"],
            metric=payload["metric"],
        ).validate()


@dataclass(frozen=True)
class FrameInvarianceObservationV231:
    observable: str
    unit: str
    value_shape: tuple[int, ...]
    base_values: tuple[complex, ...]
    translated_values: tuple[complex, ...]
    rotated_values: tuple[complex, ...]
    expected_rotated_values: tuple[complex, ...]
    translation_bohr: tuple[float, float, float]
    rotation_matrix: tuple[tuple[float, float, float], ...]
    source_artifacts: tuple[str, str, str]
    tolerance: float
    metric: str = "maximum_absolute"

    def validate(self):
        _required_text_v231("frame observable", self.observable)
        if self.unit not in _UNITS_V231:
            raise ValueError(f"unsupported evidence unit {self.unit!r}.")
        if self.metric not in _METRICS_V231:
            raise ValueError(f"unsupported evidence metric {self.metric!r}.")
        _finite_nonnegative_v231("frame tolerance", self.tolerance)
        translation = np.asarray(self.translation_bohr, dtype=float)
        rotation = np.asarray(self.rotation_matrix, dtype=float)
        if translation.shape != (3,) or not np.all(np.isfinite(translation)):
            raise ValueError("translation_bohr must be a finite three-vector.")
        if np.linalg.norm(translation) == 0.0:
            raise ValueError("frame evidence requires a nonzero translation.")
        if rotation.shape != (3, 3) or not np.all(np.isfinite(rotation)):
            raise ValueError("rotation_matrix must be finite with shape (3,3).")
        if np.linalg.norm(rotation.T @ rotation - np.eye(3), ord="fro") > 1.0e-12:
            raise ValueError("rotation_matrix must be orthogonal.")
        if abs(float(np.linalg.det(rotation)) - 1.0) > 1.0e-12:
            raise ValueError("rotation_matrix must be a proper rotation.")
        if len(self.source_artifacts) != 3 or len(set(self.source_artifacts)) != 3:
            raise ValueError("frame evidence requires three distinct source artifacts.")
        for artifact in self.source_artifacts:
            _required_text_v231("frame source artifact", artifact)
        for name, value in (
            ("base_values", self.base_values),
            ("translated_values", self.translated_values),
            ("rotated_values", self.rotated_values),
            ("expected_rotated_values", self.expected_rotated_values),
        ):
            _complex_tuple_v231(name, value, self.value_shape)
        return self

    @property
    def translation_residual(self):
        self.validate()
        return _metric_v231(
            _complex_tuple_v231("translated_values", self.translated_values, self.value_shape),
            _complex_tuple_v231("base_values", self.base_values, self.value_shape),
            self.metric,
        )

    @property
    def rotation_residual(self):
        self.validate()
        return _metric_v231(
            _complex_tuple_v231("rotated_values", self.rotated_values, self.value_shape),
            _complex_tuple_v231(
                "expected_rotated_values", self.expected_rotated_values, self.value_shape
            ),
            self.metric,
        )

    @property
    def passed(self):
        return bool(
            max(self.translation_residual, self.rotation_residual)
            <= float(self.tolerance)
        )

    def as_dict(self):
        self.validate()
        return {
            "observable": self.observable,
            "unit": self.unit,
            "value_shape": list(self.value_shape),
            "base_values": _complex_list_v231(self.base_values),
            "translated_values": _complex_list_v231(self.translated_values),
            "rotated_values": _complex_list_v231(self.rotated_values),
            "expected_rotated_values": _complex_list_v231(
                self.expected_rotated_values
            ),
            "translation_bohr": [float(value) for value in self.translation_bohr],
            "rotation_matrix": [
                [float(value) for value in row] for row in self.rotation_matrix
            ],
            "source_artifacts": list(self.source_artifacts),
            "tolerance": float(self.tolerance),
            "metric": self.metric,
            "derived_translation_residual": self.translation_residual,
            "derived_rotation_residual": self.rotation_residual,
            "passed": self.passed,
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            observable=payload["observable"],
            unit=payload["unit"],
            value_shape=tuple(payload["value_shape"]),
            base_values=_complex_from_list_v231(payload["base_values"]),
            translated_values=_complex_from_list_v231(payload["translated_values"]),
            rotated_values=_complex_from_list_v231(payload["rotated_values"]),
            expected_rotated_values=_complex_from_list_v231(
                payload["expected_rotated_values"]
            ),
            translation_bohr=tuple(payload["translation_bohr"]),
            rotation_matrix=tuple(tuple(row) for row in payload["rotation_matrix"]),
            source_artifacts=tuple(payload["source_artifacts"]),
            tolerance=payload["tolerance"],
            metric=payload["metric"],
        ).validate()


@dataclass(frozen=True)
class TrackingSpecificationV231:
    manifold_labels: tuple[str, ...]
    manifold_state_indices: tuple[tuple[int, ...], ...]
    record_edges: tuple[tuple[int, int], ...]
    overlap_threshold: float
    margin_threshold: float

    def validate(self, *, nrecord=None, nstate=None):
        labels = tuple(str(value).strip() for value in self.manifold_labels)
        if not labels or any(not value for value in labels):
            raise ValueError("tracking requires named physical manifolds.")
        if len(labels) != len(self.manifold_state_indices):
            raise ValueError("tracking labels and state-index groups disagree.")
        if len(set(labels)) != len(labels):
            raise ValueError("tracking manifold labels must be unique.")
        groups = tuple(tuple(int(value) for value in group) for group in self.manifold_state_indices)
        if any(not group or len(set(group)) != len(group) for group in groups):
            raise ValueError("tracking manifolds must contain unique state indices.")
        flattened = tuple(value for group in groups for value in group)
        if len(set(flattened)) != len(flattened) or any(value < 0 for value in flattened):
            raise ValueError("tracking manifolds must be disjoint and nonnegative.")
        if nstate is not None and set(flattened) != set(range(int(nstate))):
            raise ValueError("tracking manifolds must partition every electronic state.")
        edges = tuple((int(left), int(right)) for left, right in self.record_edges)
        if not edges or any(left == right for left, right in edges):
            raise ValueError("tracking requires non-self record edges.")
        undirected = {tuple(sorted(edge)) for edge in edges}
        if len(undirected) != len(edges):
            raise ValueError("tracking record edges must be unique up to direction.")
        if nrecord is not None:
            nrecord = int(nrecord)
            if any(
                left < 0 or right < 0 or left >= nrecord or right >= nrecord
                for left, right in edges
            ):
                raise ValueError("tracking edge contains an out-of-range record index.")
            visited = {0}
            while True:
                expanded = set(visited)
                for left, right in undirected:
                    if left in visited:
                        expanded.add(right)
                    if right in visited:
                        expanded.add(left)
                if expanded == visited:
                    break
                visited = expanded
            if visited != set(range(nrecord)):
                raise ValueError("tracking record graph must be connected.")
        _probability_v231("tracking overlap threshold", self.overlap_threshold)
        _probability_v231("tracking margin threshold", self.margin_threshold)
        return self

    def derive(self, overlaps):
        overlaps = np.asarray(overlaps, dtype=complex)
        if overlaps.ndim != 4 or overlaps.shape[0] != overlaps.shape[1]:
            raise ValueError("tracking overlap table must have shape (R,R,s,s).")
        if overlaps.shape[2] != overlaps.shape[3]:
            raise ValueError("tracking electronic overlap blocks must be square.")
        nrecord, _, nstate, _ = overlaps.shape
        self.validate(nrecord=nrecord, nstate=nstate)
        minimum_overlap = 1.0
        minimum_margin = 1.0
        rows = []
        groups = [np.asarray(group, dtype=int) for group in self.manifold_state_indices]
        for left_record, right_record in self.record_edges:
            matrix = overlaps[left_record, right_record]
            for group_index, left_group in enumerate(groups):
                right_group = groups[group_index]
                assigned = matrix[np.ix_(left_group, right_group)]
                singular_values = np.linalg.svd(assigned, compute_uv=False)
                assigned_score = float(np.min(singular_values))
                competitors = []
                for other_index, other_group in enumerate(groups):
                    if other_index == group_index:
                        continue
                    block = matrix[np.ix_(left_group, other_group)]
                    competitors.append(float(np.linalg.norm(block, ord=2)))
                competitor_score = max(competitors, default=0.0)
                margin = assigned_score - competitor_score
                minimum_overlap = min(minimum_overlap, assigned_score)
                minimum_margin = min(minimum_margin, margin)
                rows.append(
                    {
                        "left_record": int(left_record),
                        "right_record": int(right_record),
                        "manifold": self.manifold_labels[group_index],
                        "minimum_assigned_singular_value": assigned_score,
                        "maximum_competitor_spectral_norm": competitor_score,
                        "assignment_margin": margin,
                    }
                )
        return {
            "minimum_overlap": float(minimum_overlap),
            "minimum_margin": float(minimum_margin),
            "passed": bool(
                minimum_overlap >= float(self.overlap_threshold)
                and minimum_margin >= float(self.margin_threshold)
            ),
            "rows": rows,
        }

    def as_dict(self):
        self.validate()
        return {
            "manifold_labels": list(self.manifold_labels),
            "manifold_state_indices": [list(group) for group in self.manifold_state_indices],
            "record_edges": [list(edge) for edge in self.record_edges],
            "overlap_threshold": float(self.overlap_threshold),
            "margin_threshold": float(self.margin_threshold),
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            manifold_labels=tuple(payload["manifold_labels"]),
            manifold_state_indices=tuple(
                tuple(group) for group in payload["manifold_state_indices"]
            ),
            record_edges=tuple(tuple(edge) for edge in payload["record_edges"]),
            overlap_threshold=payload["overlap_threshold"],
            margin_threshold=payload["margin_threshold"],
        ).validate()


@dataclass(frozen=True)
class DerivedEvidenceBundleV231:
    reference: IndependentReferenceObservationV231
    basis: ConvergenceLadderObservationV231
    method: ConvergenceLadderObservationV231
    frame: FrameInvarianceObservationV231
    tracking: TrackingSpecificationV231

    def validate(self):
        self.reference.validate()
        self.basis.validate()
        self.method.validate()
        self.frame.validate()
        self.tracking.validate()
        if self.basis.kind != "basis" or self.method.kind != "method":
            raise ValueError("derived evidence requires basis and method ladders.")
        return self

    def derive_v230(self, overlaps):
        self.validate()
        tracking = self.tracking.derive(overlaps)
        return MolecularSOCValidationEvidenceV230(
            independent_reference_id=self.reference.reference_id,
            independent_reference_error=self.reference.error,
            independent_reference_tolerance=float(self.reference.tolerance),
            basis_levels=tuple(self.basis.labels),
            basis_changes=tuple(self.basis.changes),
            basis_tolerance=float(self.basis.tolerance),
            method_levels=tuple(self.method.labels),
            method_changes=tuple(self.method.changes),
            method_tolerance=float(self.method.tolerance),
            translation_residual=self.frame.translation_residual,
            rotation_residual=self.frame.rotation_residual,
            frame_invariance_tolerance=float(self.frame.tolerance),
            tracking_minimum_overlap=tracking["minimum_overlap"],
            tracking_minimum_margin=tracking["minimum_margin"],
            tracking_overlap_threshold=float(self.tracking.overlap_threshold),
            tracking_margin_threshold=float(self.tracking.margin_threshold),
        ).validate()

    def as_dict(self):
        self.validate()
        return {
            "reference": self.reference.as_dict(),
            "basis": self.basis.as_dict(),
            "method": self.method.as_dict(),
            "frame": self.frame.as_dict(),
            "tracking": self.tracking.as_dict(),
        }

    @classmethod
    def from_dict(cls, payload):
        return cls(
            reference=IndependentReferenceObservationV231.from_dict(payload["reference"]),
            basis=ConvergenceLadderObservationV231.from_dict(payload["basis"]),
            method=ConvergenceLadderObservationV231.from_dict(payload["method"]),
            frame=FrameInvarianceObservationV231.from_dict(payload["frame"]),
            tracking=TrackingSpecificationV231.from_dict(payload["tracking"]),
        ).validate()
