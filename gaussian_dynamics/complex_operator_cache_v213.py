"""Fingerprint-safe disk cache for fixed-frame complex electronic operators."""

from pathlib import Path
import hashlib
import json
import numpy as np

from .electronic_contract_v213 import (
    ElectronicOperatorProvenanceV213,
    validate_electronic_contract_v213,
)
from .electronic_operator_v21 import (
    ElectronicOperatorPointV21,
    ElectronicOperatorSnapshotV21,
)
from .matrix_invariants_v213 import require_residual_v213, scaled_matrix_residual_v213


def _json_encode(value):
    if isinstance(value, np.generic):
        return _json_encode(value.item())
    if isinstance(value, complex):
        return {"__complex__": [float(value.real), float(value.imag)]}
    if isinstance(value, np.ndarray):
        return _json_encode(value.tolist())
    if isinstance(value, dict):
        return {str(key): _json_encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"cache metadata contains unsupported {type(value).__name__}.")


def _json_decode(value):
    if isinstance(value, dict):
        if set(value) == {"__complex__"}:
            real, imag = value["__complex__"]
            return complex(real, imag)
        return {key: _json_decode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_decode(item) for item in value]
    return value


class FixedFrameComplexOperatorCacheV213:
    """Cache H, K, D, mass, and state vectors with a mandatory model fingerprint.

    This cache deliberately does not accept moving-frame wavefunction snapshots.  Those
    require a separate ab-initio overlap/restart design and must not be confused with a
    fixed spin-diabatic or fixed general electronic frame.
    """

    format_version = "v0.21.3-fixed-frame-complex-operator-1"

    def __init__(
        self,
        provider,
        directory,
        provenance: ElectronicOperatorProvenanceV213,
        *,
        namespace="default",
        coordinate_digits=12,
    ):
        self.provider = provider
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.provenance = provenance.validate()
        if not self.provenance.model_space.fixed_frame:
            raise ValueError("FixedFrameComplexOperatorCacheV213 requires a fixed frame.")
        self.namespace = str(namespace)
        if not self.namespace:
            raise ValueError("cache namespace cannot be empty.")
        self.coordinate_digits = int(coordinate_digits)
        if not 0 <= self.coordinate_digits <= 15:
            raise ValueError("coordinate_digits must lie between 0 and 15.")
        self.hits = 0
        self.misses = 0

    @property
    def provider_fingerprint(self):
        return self.provenance.fingerprint()

    def _key(self, q):
        payload = {
            "format": self.format_version,
            "namespace": self.namespace,
            "provider_fingerprint": self.provider_fingerprint,
            "coordinate_digits": self.coordinate_digits,
            "q": np.asarray(q, dtype=float).round(self.coordinate_digits).tolist(),
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(raw).hexdigest()

    def _paths(self, q):
        key = self._key(q)
        return self.directory / f"{key}.npz", self.directory / f"{key}.json"

    def _point_with_contract(self, point, *, cache_status):
        metadata = {
            **dict(point.metadata),
            "v213_electronic_contract": self.provenance.as_dict(),
            "v213_provenance_fingerprint": self.provider_fingerprint,
            "v213_cache": {
                "format": self.format_version,
                "namespace": self.namespace,
                "status": cache_status,
            },
        }
        contracted = ElectronicOperatorPointV21(
            q=np.asarray(point.q, dtype=float).copy(),
            H=np.asarray(point.H, dtype=complex).copy(),
            dH_dq=np.asarray(point.dH_dq, dtype=complex).copy(),
            connection_q=np.asarray(point.connection_q, dtype=complex).copy(),
            mass_matrix_q_au=np.asarray(point.mass_matrix_q_au, dtype=float).copy(),
            metadata=metadata,
        )
        return validate_electronic_contract_v213(contracted, self.provenance)

    def _load(self, npz_path, json_path, requested_q):
        with np.load(npz_path, allow_pickle=False) as arrays:
            metadata_payload = json.loads(json_path.read_text(encoding="utf-8"))
            if metadata_payload["format"] != self.format_version:
                raise ValueError("complex-operator cache format mismatch.")
            if metadata_payload["provider_fingerprint"] != self.provider_fingerprint:
                raise ValueError("complex-operator cache provenance mismatch.")
            if metadata_payload["namespace"] != self.namespace:
                raise ValueError("complex-operator cache namespace mismatch.")
            if int(metadata_payload["coordinate_digits"]) != self.coordinate_digits:
                raise ValueError("complex-operator cache coordinate policy mismatch.")
            cached_q = arrays["q"]
            if (
                cached_q.shape != requested_q.shape
                or not np.array_equal(
                    cached_q.round(self.coordinate_digits),
                    requested_q.round(self.coordinate_digits),
                )
            ):
                raise ValueError("complex-operator cache coordinate mismatch.")
            point = ElectronicOperatorPointV21(
                q=cached_q,
                H=arrays["H"],
                dH_dq=arrays["dH_dq"],
                connection_q=arrays["connection_q"],
                mass_matrix_q_au=arrays["mass_matrix_q_au"],
                metadata=_json_decode(metadata_payload["point_metadata"]),
            )
            vectors = arrays["state_vectors"] if "state_vectors" in arrays.files else None
        point = self._point_with_contract(point, cache_status="hit")
        return ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=vectors,
            metadata={"cache_status": "hit", "provider_fingerprint": self.provider_fingerprint},
        ).validate()

    def _store(self, npz_path, json_path, snapshot):
        temp_npz = npz_path.with_suffix(".npz.tmp")
        temp_json = json_path.with_suffix(".json.tmp")
        arrays = {
            "q": snapshot.point.q,
            "H": snapshot.point.H,
            "dH_dq": snapshot.point.dH_dq,
            "connection_q": snapshot.point.connection_q,
            "mass_matrix_q_au": snapshot.point.mass_matrix_q_au,
        }
        if snapshot.state_vectors is not None:
            arrays["state_vectors"] = snapshot.state_vectors
        with temp_npz.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
        payload = {
            "format": self.format_version,
            "namespace": self.namespace,
            "provider_fingerprint": self.provider_fingerprint,
            "coordinate_digits": self.coordinate_digits,
            "point_metadata": _json_encode(snapshot.point.metadata),
        }
        temp_json.write_text(
            json.dumps(payload, sort_keys=True, indent=2), encoding="utf-8"
        )
        temp_npz.replace(npz_path)
        temp_json.replace(json_path)

    def evaluate_snapshot(self, q):
        raw_q = np.asarray(q)
        if np.iscomplexobj(raw_q) and np.any(np.imag(raw_q) != 0.0):
            raise ValueError("cache coordinates must be real.")
        q = np.asarray(np.real(raw_q), dtype=float)
        if q.ndim != 1 or not np.all(np.isfinite(q)):
            raise ValueError("cache coordinates must be a finite vector.")
        npz_path, json_path = self._paths(q)
        if npz_path.exists() != json_path.exists():
            raise RuntimeError("complex-operator cache entry is incomplete.")
        if npz_path.exists():
            self.hits += 1
            return self._load(npz_path, json_path, q)

        self.misses += 1
        raw = self.provider.evaluate_snapshot(q)
        if raw.wavefunction_snapshot is not None:
            raise ValueError(
                "moving-frame wavefunction snapshots need a dedicated molecular cache."
            )
        point = self._point_with_contract(raw.point, cache_status="miss")
        vectors = (
            np.eye(point.nstate, dtype=complex)
            if raw.state_vectors is None
            else np.asarray(raw.state_vectors, dtype=complex).copy()
        )
        snapshot = ElectronicOperatorSnapshotV21(
            point=point,
            state_vectors=vectors,
            metadata={"cache_status": "miss", "provider_fingerprint": self.provider_fingerprint},
        ).validate()
        self._store(npz_path, json_path, snapshot)
        return snapshot

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    @staticmethod
    def snapshot_overlap(left, right):
        overlap = left.state_vectors.conj().T @ right.state_vectors
        require_residual_v213(
            "fixed-frame cross-geometry overlap identity",
            scaled_matrix_residual_v213(
                overlap, np.eye(left.point.nstate, dtype=complex)
            ),
            1.0e-10,
        )
        return overlap

    def diagnostics_dict(self):
        return {
            "format": self.format_version,
            "namespace": self.namespace,
            "provider_fingerprint": self.provider_fingerprint,
            "coordinate_digits": int(self.coordinate_digits),
            "hits": int(self.hits),
            "misses": int(self.misses),
        }
