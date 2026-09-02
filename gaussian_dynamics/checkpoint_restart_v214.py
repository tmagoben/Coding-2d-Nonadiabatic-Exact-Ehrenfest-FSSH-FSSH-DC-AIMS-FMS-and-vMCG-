"""Provenance-checked deterministic checkpoint/restart for v0.21.4.

The checkpoint contains the full Gaussian block state, density-guide state, and
sparse-graph hysteresis state.  It deliberately stores no provider wavefunction
objects: those are reconstructed at the checkpoint geometries through the declared
provider and verified with the electronic provenance fingerprint.
"""

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import numpy as np

from .block_sparse_molecular_v21 import BlockMolecularTBFV21
from .density_guidance_v213 import (
    BlockDensityMatrixGuidanceV213,
    validate_guide_density_v213,
)
from .electronic_contract_v213 import validate_electronic_contract_v213
from .self_consistent_block_v212 import run_self_consistent_block_dynamics_v212
from .self_consistent_block_v213 import SelfConsistentBlockSettingsV213


CHECKPOINT_FORMAT_V214 = "v0.21.4-self-consistent-block-checkpoint-1"
_ARRAY_NAMES_V214 = (
    "uids",
    "q",
    "p",
    "A",
    "coefficients",
    "guide_mask",
    "guide_densities",
    "guidance_counters",
    "active_uid_edges",
)


def _canonical_value_v214(value):
    if isinstance(value, np.generic):
        return _canonical_value_v214(value.item())
    if isinstance(value, np.ndarray):
        return _canonical_value_v214(value.tolist())
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise TypeError("checkpoint metadata dictionary keys must be strings.")
        return {key: _canonical_value_v214(value[key]) for key in sorted(value)}
    if isinstance(value, (tuple, list)):
        return [_canonical_value_v214(item) for item in value]
    if isinstance(value, complex):
        if not np.isfinite(value.real) or not np.isfinite(value.imag):
            raise ValueError("checkpoint metadata cannot contain non-finite values.")
        return {"__complex__": [float(value.real), float(value.imag)]}
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("checkpoint metadata cannot contain non-finite values.")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported checkpoint metadata type {type(value).__name__}.")


def _canonical_json_v214(value):
    return json.dumps(
        _canonical_value_v214(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def settings_fingerprint_v214(settings):
    settings = settings.validate()
    payload = {
        "release": "v0.21.4",
        "control": asdict(settings),
        "semantics": {
            "guidance": "transported density matrix",
            "trial_state": "rollback then commit accepted endpoint",
            "graph_state": "uid-edge hysteresis restored",
        },
    }
    return hashlib.sha256(_canonical_json_v214(payload).encode("utf-8")).hexdigest()


def _integrity_digest_v214(manifest_without_digest, arrays):
    digest = hashlib.sha256()
    digest.update(_canonical_json_v214(manifest_without_digest).encode("utf-8"))
    for name in _ARRAY_NAMES_V214:
        array = np.ascontiguousarray(arrays[name])
        descriptor = {
            "name": name,
            "dtype": array.dtype.str,
            "shape": list(array.shape),
        }
        digest.update(_canonical_json_v214(descriptor).encode("utf-8"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


@dataclass(frozen=True)
class SelfConsistentBlockSettingsV214(SelfConsistentBlockSettingsV213):
    checkpoint_identity_policy: str = "strict"

    def validate(self):
        super().validate()
        if self.checkpoint_identity_policy != "strict":
            raise ValueError("v0.21.4 supports only strict checkpoint identity.")
        return self


@dataclass(frozen=True)
class SelfConsistentBlockCheckpointV214:
    format_version: str
    step: int
    time: float
    dt: float
    nstate: int
    nuclear_dimension: int
    provider_fingerprint: str
    settings_fingerprint: str
    uids: np.ndarray
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray
    coefficients: np.ndarray
    guide_mask: np.ndarray
    guide_densities: np.ndarray
    guidance_counters: np.ndarray
    active_uid_edges: np.ndarray
    integrity_digest: str

    def _manifest_without_digest(self):
        return {
            "format_version": self.format_version,
            "step": int(self.step),
            "time": float(self.time),
            "dt": float(self.dt),
            "nstate": int(self.nstate),
            "nuclear_dimension": int(self.nuclear_dimension),
            "provider_fingerprint": self.provider_fingerprint,
            "settings_fingerprint": self.settings_fingerprint,
        }

    def arrays_dict(self):
        return {
            "uids": np.asarray(self.uids),
            "q": np.asarray(self.q),
            "p": np.asarray(self.p),
            "A": np.asarray(self.A),
            "coefficients": np.asarray(self.coefficients),
            "guide_mask": np.asarray(self.guide_mask),
            "guide_densities": np.asarray(self.guide_densities),
            "guidance_counters": np.asarray(self.guidance_counters),
            "active_uid_edges": np.asarray(self.active_uid_edges),
        }

    def computed_integrity_digest(self):
        return _integrity_digest_v214(
            self._manifest_without_digest(), self.arrays_dict()
        )

    def validate(
        self,
        *,
        expected_provider_fingerprint=None,
        expected_settings_fingerprint=None,
    ):
        if self.format_version != CHECKPOINT_FORMAT_V214:
            raise ValueError("unsupported v0.21.4 checkpoint format.")
        if int(self.step) != self.step or int(self.step) < 0:
            raise ValueError("checkpoint step must be a nonnegative integer.")
        if not np.isfinite(self.dt) or float(self.dt) <= 0.0:
            raise ValueError("checkpoint dt must be finite and positive.")
        if not np.isfinite(self.time) or float(self.time) < 0.0:
            raise ValueError("checkpoint time must be finite and nonnegative.")
        expected_time = int(self.step) * float(self.dt)
        if abs(float(self.time) - expected_time) > 1.0e-12 * max(1.0, expected_time):
            raise ValueError("checkpoint time is inconsistent with step*dt.")
        if int(self.nstate) != self.nstate or int(self.nstate) < 1:
            raise ValueError("checkpoint nstate must be a positive integer.")
        if (
            int(self.nuclear_dimension) != self.nuclear_dimension
            or int(self.nuclear_dimension) < 1
        ):
            raise ValueError("checkpoint nuclear dimension must be positive.")
        for name, fingerprint in (
            ("provider", self.provider_fingerprint),
            ("settings", self.settings_fingerprint),
            ("integrity", self.integrity_digest),
        ):
            if (
                not isinstance(fingerprint, str)
                or len(fingerprint) != 64
                or any(character not in "0123456789abcdef" for character in fingerprint)
            ):
                raise ValueError(f"checkpoint {name} fingerprint is not SHA-256 hex.")

        uids = np.asarray(self.uids)
        q = np.asarray(self.q, dtype=float)
        p = np.asarray(self.p, dtype=float)
        widths = np.asarray(self.A, dtype=float)
        coefficients = np.asarray(self.coefficients, dtype=complex)
        mask = np.asarray(self.guide_mask)
        densities = np.asarray(self.guide_densities, dtype=complex)
        counters = np.asarray(self.guidance_counters)
        edges = np.asarray(self.active_uid_edges)
        nstate = int(self.nstate)
        dimension = int(self.nuclear_dimension)
        if uids.ndim != 1 or len(uids) < 1:
            raise ValueError("checkpoint must contain at least one Gaussian uid.")
        if not np.issubdtype(uids.dtype, np.integer):
            raise ValueError("checkpoint Gaussian uids must have integer dtype.")
        uids = uids.astype(np.int64, copy=False)
        if len(set(uids.tolist())) != len(uids):
            raise ValueError("checkpoint Gaussian uids must be unique.")
        n = len(uids)
        if q.shape != (n, dimension) or p.shape != (n, dimension):
            raise ValueError("checkpoint q/p arrays have incompatible shapes.")
        if widths.shape != (n, dimension, dimension):
            raise ValueError("checkpoint width matrices have incompatible shape.")
        if coefficients.shape != (n * nstate,):
            raise ValueError("checkpoint coefficient vector has incompatible shape.")
        if mask.shape != (n,) or mask.dtype != np.dtype(bool):
            raise ValueError("checkpoint guide_mask must be a Boolean vector.")
        if densities.shape != (n, nstate, nstate):
            raise ValueError("checkpoint guide densities have incompatible shape.")
        if counters.shape != (6,) or not np.issubdtype(counters.dtype, np.integer):
            raise ValueError("checkpoint guidance counters must be six integers.")
        if np.any(counters < 0):
            raise ValueError("checkpoint guidance counters cannot be negative.")
        if edges.ndim != 2 or edges.shape[1:] != (2,):
            raise ValueError("checkpoint active_uid_edges must have shape (E,2).")
        if not np.issubdtype(edges.dtype, np.integer):
            raise ValueError("checkpoint active graph uids must have integer dtype.")
        if not all(
            np.all(np.isfinite(array))
            for array in (q, p, widths, coefficients, densities)
        ):
            raise ValueError("checkpoint numerical arrays contain non-finite data.")
        for index, width in enumerate(widths):
            if not np.allclose(width, width.T, rtol=0.0, atol=1.0e-12):
                raise ValueError(f"checkpoint width matrix {index} is not symmetric.")
            if float(np.min(np.linalg.eigvalsh(width))) <= 0.0:
                raise ValueError(f"checkpoint width matrix {index} is not positive definite.")
        for index in np.flatnonzero(mask):
            validate_guide_density_v213(densities[index], nstate)
        if np.linalg.norm(densities[~mask]) != 0.0:
            raise ValueError("untracked checkpoint guide densities must be exact zeros.")
        live = set(uids.tolist())
        canonical_edges = []
        for raw_a, raw_b in edges.astype(np.int64, copy=False):
            if raw_a == raw_b or raw_a not in live or raw_b not in live:
                raise ValueError("checkpoint graph edge references invalid Gaussian uids.")
            edge = (int(raw_a), int(raw_b))
            canonical = (min(edge), max(edge))
            if edge != canonical:
                raise ValueError("checkpoint graph edges must use ascending uid order.")
            canonical_edges.append(canonical)
        if len(set(canonical_edges)) != len(canonical_edges):
            raise ValueError("checkpoint graph edges must be unique.")
        if canonical_edges != sorted(canonical_edges):
            raise ValueError("checkpoint graph edges must be canonical and sorted.")
        if expected_provider_fingerprint is not None:
            if self.provider_fingerprint != expected_provider_fingerprint:
                raise ValueError("checkpoint provider provenance fingerprint mismatch.")
        if expected_settings_fingerprint is not None:
            if self.settings_fingerprint != expected_settings_fingerprint:
                raise ValueError("checkpoint propagation-settings fingerprint mismatch.")
        if self.computed_integrity_digest() != self.integrity_digest:
            raise ValueError("checkpoint integrity digest mismatch.")
        return self

    def basis(self):
        self.validate()
        return tuple(
            BlockMolecularTBFV21(
                int(self.uids[index]),
                np.asarray(self.q[index], dtype=float).copy(),
                np.asarray(self.p[index], dtype=float).copy(),
                np.asarray(self.A[index], dtype=float).copy(),
            )
            for index in range(len(self.uids))
        )

    @classmethod
    def create(
        cls,
        *,
        step,
        dt,
        provider_fingerprint,
        settings_fingerprint,
        basis,
        coefficients,
        nstate,
        guidance,
        active_uid_edges=(),
    ):
        basis = tuple(basis)
        if not basis:
            raise ValueError("cannot checkpoint an empty Gaussian basis.")
        nstate = int(nstate)
        uids = np.asarray([item.uid for item in basis], dtype=np.int64)
        q = np.asarray([item.q for item in basis], dtype=float)
        p = np.asarray([item.p for item in basis], dtype=float)
        widths = np.asarray([item.A for item in basis], dtype=float)
        coefficients = np.asarray(coefficients, dtype=complex).copy()
        guidance_state = guidance.checkpoint_state()
        guide_mask = np.asarray(
            [int(uid) in guidance_state["densities"] for uid in uids], dtype=bool
        )
        guide_densities = np.zeros((len(basis), nstate, nstate), dtype=complex)
        for index, uid in enumerate(uids):
            if guide_mask[index]:
                guide_densities[index] = guidance_state["densities"][int(uid)]
        counters = np.asarray(guidance_state["counters"], dtype=np.int64)
        edges = np.asarray(tuple(active_uid_edges), dtype=np.int64)
        if edges.size == 0:
            edges = np.empty((0, 2), dtype=np.int64)
        manifest = {
            "format_version": CHECKPOINT_FORMAT_V214,
            "step": int(step),
            "time": float(step) * float(dt),
            "dt": float(dt),
            "nstate": int(nstate),
            "nuclear_dimension": int(q.shape[1]),
            "provider_fingerprint": str(provider_fingerprint),
            "settings_fingerprint": str(settings_fingerprint),
        }
        arrays = {
            "uids": uids,
            "q": q,
            "p": p,
            "A": widths,
            "coefficients": coefficients,
            "guide_mask": guide_mask,
            "guide_densities": guide_densities,
            "guidance_counters": counters,
            "active_uid_edges": edges,
        }
        checkpoint = cls(
            **manifest,
            **arrays,
            integrity_digest=_integrity_digest_v214(manifest, arrays),
        )
        return checkpoint.validate()


def save_self_consistent_checkpoint_v214(path, checkpoint):
    checkpoint = checkpoint.validate()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    manifest = {
        **checkpoint._manifest_without_digest(),
        "integrity_digest": checkpoint.integrity_digest,
    }
    with temporary.open("wb") as handle:
        np.savez_compressed(
            handle,
            manifest_json=np.asarray(_canonical_json_v214(manifest)),
            **checkpoint.arrays_dict(),
        )
    temporary.replace(path)
    return path


def load_self_consistent_checkpoint_v214(
    path,
    *,
    expected_provider_fingerprint=None,
    expected_settings_fingerprint=None,
):
    path = Path(path)
    with np.load(path, allow_pickle=False) as archive:
        expected = {"manifest_json", *_ARRAY_NAMES_V214}
        if set(archive.files) != expected:
            raise ValueError("checkpoint archive has missing or unexpected arrays.")
        manifest_array = archive["manifest_json"]
        if manifest_array.shape != () or manifest_array.dtype.kind not in {"U", "S"}:
            raise ValueError("checkpoint manifest must be one JSON string.")
        manifest_text = manifest_array.item()
        if isinstance(manifest_text, bytes):
            manifest_text = manifest_text.decode("utf-8")
        manifest = json.loads(manifest_text)
        required_manifest = {
            "format_version",
            "step",
            "time",
            "dt",
            "nstate",
            "nuclear_dimension",
            "provider_fingerprint",
            "settings_fingerprint",
            "integrity_digest",
        }
        if set(manifest) != required_manifest:
            raise ValueError("checkpoint manifest fields do not match the format.")
        arrays = {name: archive[name].copy() for name in _ARRAY_NAMES_V214}
    checkpoint = SelfConsistentBlockCheckpointV214(
        **manifest,
        **arrays,
    )
    return checkpoint.validate(
        expected_provider_fingerprint=expected_provider_fingerprint,
        expected_settings_fingerprint=expected_settings_fingerprint,
    )


def _restore_guidance_v214(checkpoint, provider, settings):
    guidance = BlockDensityMatrixGuidanceV213(settings.guidance)
    basis = checkpoint.basis()
    for index, item in enumerate(basis):
        guide_density = (
            checkpoint.guide_densities[index]
            if checkpoint.guide_mask[index]
            else None
        )
        guidance.on_insert(item, provider, guide_density=guide_density)
    state = guidance.checkpoint_state()
    state["counters"] = tuple(int(value) for value in checkpoint.guidance_counters)
    guidance.restore_state(state)
    return guidance


def run_self_consistent_block_dynamics_v214(
    provider,
    provenance,
    *,
    initial_basis=None,
    C0=None,
    checkpoint=None,
    dt=None,
    steps=20,
    settings=SelfConsistentBlockSettingsV214(),
    store_every=5,
    adaptation_policy=None,
    provider_numerical_fingerprint=None,
):
    """Run a fresh or resumed density-guided trajectory and return a new checkpoint."""
    settings = settings.validate()
    provenance = provenance.validate()
    provenance_fingerprint = provenance.fingerprint()
    if provider_numerical_fingerprint is None:
        provider_fingerprint = provenance_fingerprint
    else:
        provider_fingerprint = str(provider_numerical_fingerprint)
        if len(provider_fingerprint) != 64 or any(
            character not in "0123456789abcdef"
            for character in provider_fingerprint
        ):
            raise ValueError(
                "provider_numerical_fingerprint must be a lowercase SHA-256 digest."
            )
    settings_fingerprint = settings_fingerprint_v214(settings)
    if checkpoint is None:
        if initial_basis is None or C0 is None:
            raise ValueError("a fresh v0.21.4 run requires initial_basis and C0.")
        basis = tuple(item.copy() for item in initial_basis)
        coefficients = np.asarray(C0, dtype=complex).copy()
        offset_step = 0
        offset_time = 0.0
        dt = 0.002 if dt is None else float(dt)
        guidance = BlockDensityMatrixGuidanceV213(settings.guidance)
        graph_edges = None
    else:
        if initial_basis is not None or C0 is not None:
            raise ValueError("a resumed run cannot also specify initial_basis or C0.")
        checkpoint = checkpoint.validate(
            expected_provider_fingerprint=provider_fingerprint,
            expected_settings_fingerprint=settings_fingerprint,
        )
        if dt is not None and float(dt) != float(checkpoint.dt):
            raise ValueError("restart dt differs from the checkpoint dt.")
        dt = float(checkpoint.dt)
        basis = checkpoint.basis()
        coefficients = np.asarray(checkpoint.coefficients, dtype=complex).copy()
        offset_step = int(checkpoint.step)
        offset_time = float(checkpoint.time)
        guidance = _restore_guidance_v214(checkpoint, provider, settings)
        graph_edges = tuple(map(tuple, checkpoint.active_uid_edges.tolist()))

    if not basis:
        raise ValueError("v0.21.4 propagation requires a nonempty Gaussian basis.")
    if not np.isfinite(dt) or float(dt) <= 0.0:
        raise ValueError("v0.21.4 propagation dt must be finite and positive.")
    nstate = provenance.model_space.nstate
    for item in basis:
        snapshot = provider.evaluate_snapshot(item.q).validate()
        validate_electronic_contract_v213(snapshot.point, provenance)
        if snapshot.point.nstate != nstate:
            raise ValueError("provider dimension differs from checkpoint/model space.")
        emitted_fingerprint = snapshot.point.metadata.get(
            "v213_provenance_fingerprint"
        )
        if emitted_fingerprint != provenance_fingerprint:
            raise ValueError(
                "v0.21.4 propagation requires the provider to emit the declared "
                "provenance fingerprint."
            )

    shifted_policy = adaptation_policy
    if adaptation_policy is not None and offset_step:
        def shifted_policy(local_step, local_basis, local_C, local_S):
            return adaptation_policy(
                offset_step + int(local_step), local_basis, local_C, local_S
            )

    output = run_self_consistent_block_dynamics_v212(
        basis,
        coefficients,
        provider,
        dt=dt,
        steps=steps,
        settings=settings.legacy_control_settings(),
        store_every=store_every,
        adaptation_policy=shifted_policy,
        guidance_engine=guidance,
        graph_active_uid_edges=graph_edges,
    )
    for record in output["records"]:
        record["step"] += offset_step
        record["time"] += offset_time
    for row in output["corrector_history"]:
        row["step"] += offset_step
    for event in output["adaptation_events"]:
        event["step"] += offset_step
    total_step = offset_step + int(steps)
    final_checkpoint = SelfConsistentBlockCheckpointV214.create(
        step=total_step,
        dt=dt,
        provider_fingerprint=provider_fingerprint,
        settings_fingerprint=settings_fingerprint,
        basis=output["final_basis"],
        coefficients=output["final_coefficients"],
        nstate=nstate,
        guidance=guidance,
        active_uid_edges=output["final_active_uid_edges"],
    )
    output["settings"] = {
        "dt": float(dt),
        "segment_steps": int(steps),
        "initial_step": int(offset_step),
        "final_step": int(total_step),
        "control": asdict(settings),
        "settings_fingerprint": settings_fingerprint,
        "provider_fingerprint": provider_fingerprint,
    }
    output["release_path"] = "v0.21.4"
    output["restart_source"] = checkpoint is not None
    output["checkpoint"] = final_checkpoint
    return output
