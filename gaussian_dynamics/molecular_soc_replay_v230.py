"""Deterministic, integrity-checked molecular SOC snapshot replay for v0.23.0."""

from dataclasses import asdict, dataclass
from io import BytesIO
import hashlib
import json
import os
from pathlib import Path
import zipfile
import numpy as np

from .analytic_soc_models_v220 import SOCOperatorComponentsV220
from .electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
    compose_electronic_operator_v213,
)
from .electronic_operator_v21 import ElectronicOperatorSnapshotV21
from .molecular_soc_contract_v230 import (
    MolecularSOCAdmissionContractV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCCapabilitiesV230,
    MolecularSOCValidationEvidenceV230,
    provenance_with_molecular_soc_contract_v230,
)
from .soc_admission_v221 import SOCSymmetryContractV221


REPLAY_MANIFEST_NAME_V230 = "molecular_soc_manifest_v230.json"
REPLAY_ARRAYS_NAME_V230 = "molecular_soc_arrays_v230.npz"


def _sha256_file_v230(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_json_bytes_v230(value):
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("utf-8")


def _write_deterministic_npz_v230(path, arrays):
    """Write sorted NPY members with fixed ZIP metadata for byte reproducibility."""
    path = Path(path)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in sorted(arrays):
            buffer = BytesIO()
            np.lib.format.write_array(
                buffer,
                np.ascontiguousarray(arrays[name]),
                allow_pickle=False,
            )
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, buffer.getvalue(), compress_type=zipfile.ZIP_DEFLATED)


def _decode_canonical_complex_v230(value):
    if isinstance(value, dict):
        if set(value) == {"__complex__"}:
            pair = value["__complex__"]
            if not isinstance(pair, list) or len(pair) != 2:
                raise ValueError("invalid canonical complex value in replay manifest.")
            return complex(float(pair[0]), float(pair[1]))
        return {
            str(key): _decode_canonical_complex_v230(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_decode_canonical_complex_v230(item) for item in value]
    return value


def _provenance_from_dict_v230(payload):
    space = payload["model_space"]
    model_space = ElectronicModelSpaceV213(
        name=space["name"],
        representation=space["representation"],
        states=tuple(
            ElectronicStateDescriptorV213(
                label=state["label"],
                source_root=state.get("source_root"),
                multiplicity=state.get("multiplicity"),
                component=state.get("component"),
                charge=state.get("charge"),
            )
            for state in space["states"]
        ),
        fixed_dimension=space["fixed_dimension"],
        full_electronic_blocks=space["full_electronic_blocks"],
        complete_multiplets=space["complete_multiplets"],
        energy_unit=space["energy_unit"],
        coordinate_unit=space["coordinate_unit"],
    )
    return ElectronicOperatorProvenanceV213(
        model_name=payload["model_name"],
        model_version=payload["model_version"],
        model_space=model_space,
        spin_free_method=payload["spin_free_method"],
        soc_enabled=payload["soc_enabled"],
        soc_method=payload["soc_method"],
        scalar_relativistic_method=payload["scalar_relativistic_method"],
        derivative_method=payload["derivative_method"],
        parameters=_decode_canonical_complex_v230(payload["parameters"]),
    ).validate()


def _contract_from_dict_v230(payload):
    capabilities = MolecularSOCCapabilitiesV230(
        **{
            name: payload["capabilities"][name]
            for name in (
                "static_soc",
                "spin_free_derivatives",
                "soc_derivatives",
                "derivative_connections",
                "cross_geometry_overlaps",
                "deterministic_replay",
                "analytic_soc_derivatives",
            )
        }
    )
    identity_payload = payload["identity"]
    identity = MolecularSOCBackendIdentityV230(
        **{
            name: identity_payload[name]
            for name in (
                "backend_name",
                "backend_version",
                "source_kind",
                "electronic_method",
                "basis",
                "charge",
                "electron_count",
                "soc_operator",
                "scalar_relativistic_method",
                "derivative_method",
                "active_space",
                "molecule_name",
                "atom_symbols",
                "isotope_masses_amu",
                "reference_geometry_bohr",
                "calculation_input_sha256",
                "environment_sha256",
                "geometry_unit",
                "energy_unit",
                "derivative_unit",
                "extra",
            )
        }
    )
    evidence_payload = payload["evidence"]
    evidence = MolecularSOCValidationEvidenceV230(
        independent_reference_id=evidence_payload["independent_reference_id"],
        independent_reference_error=evidence_payload["independent_reference_error"],
        independent_reference_tolerance=evidence_payload[
            "independent_reference_tolerance"
        ],
        basis_levels=tuple(evidence_payload["basis_levels"]),
        basis_changes=tuple(evidence_payload["basis_changes"]),
        basis_tolerance=evidence_payload["basis_tolerance"],
        method_levels=tuple(evidence_payload["method_levels"]),
        method_changes=tuple(evidence_payload["method_changes"]),
        method_tolerance=evidence_payload["method_tolerance"],
        translation_residual=evidence_payload["translation_residual"],
        rotation_residual=evidence_payload["rotation_residual"],
        frame_invariance_tolerance=evidence_payload["frame_invariance_tolerance"],
        tracking_minimum_overlap=evidence_payload["tracking_minimum_overlap"],
        tracking_minimum_margin=evidence_payload["tracking_minimum_margin"],
        tracking_overlap_threshold=evidence_payload["tracking_overlap_threshold"],
        tracking_margin_threshold=evidence_payload["tracking_margin_threshold"],
    )
    return MolecularSOCAdmissionContractV230(
        capabilities=capabilities,
        identity=identity,
        evidence=evidence,
        state_tracking_policy=payload["state_tracking_policy"],
        coordinate_definition=payload["coordinate_definition"],
        all_electronic_calculations_converged=payload[
            "all_electronic_calculations_converged"
        ],
    ).validate()


@dataclass(frozen=True)
class ReplayWavefunctionTokenV230:
    dataset_fingerprint: str
    record_index: int


@dataclass(frozen=True)
class MolecularSOCReplayOverlapDiagnosticsV232:
    """Necessary Hilbert-space consistency checks for replay overlaps.

    Cross-geometry overlaps between two finite orthonormal state manifolds are
    contractions, not generally unitary matrices.  Their singular values may
    therefore be smaller than one, but never larger than one.  Self overlaps
    remain the identity and reversing a geometry pair takes the adjoint.
    """

    record_count: int
    state_count: int
    unordered_cross_geometry_pair_count: int
    maximum_self_identity_residual: float
    maximum_reciprocity_residual: float
    minimum_cross_geometry_singular_value: float | None
    maximum_cross_geometry_singular_value: float | None
    maximum_contraction_excess: float

    def as_dict(self):
        return asdict(self)


def _replay_overlap_diagnostics_v232(overlaps):
    overlaps = np.asarray(overlaps, dtype=complex)
    if overlaps.ndim != 4:
        raise ValueError("replay overlap diagnostics require a rank-four table.")
    nrecord_left, nrecord_right, nstate_left, nstate_right = overlaps.shape
    if (
        nrecord_left < 1
        or nrecord_left != nrecord_right
        or nstate_left < 1
        or nstate_left != nstate_right
    ):
        raise ValueError(
            "replay overlap diagnostics require square record and state dimensions."
        )
    if not np.all(np.isfinite(overlaps)):
        raise ValueError("replay overlaps contain non-finite data.")

    identity = np.eye(nstate_left, dtype=complex)
    self_residuals = [
        float(np.linalg.norm(overlaps[index, index] - identity, ord="fro"))
        for index in range(nrecord_left)
    ]
    reciprocity_residuals = []
    cross_geometry_singular_values = []
    for left in range(nrecord_left):
        for right in range(left + 1, nrecord_left):
            reciprocity_residuals.append(
                float(
                    np.linalg.norm(
                        overlaps[left, right]
                        - overlaps[right, left].conj().T,
                        ord="fro",
                    )
                )
            )
            # Inspect both ordered blocks so an inconsistent reverse block
            # cannot hide an expansive singular value behind failed reciprocity.
            cross_geometry_singular_values.extend(
                np.linalg.svd(overlaps[left, right], compute_uv=False).tolist()
            )
            cross_geometry_singular_values.extend(
                np.linalg.svd(overlaps[right, left], compute_uv=False).tolist()
            )

    if cross_geometry_singular_values:
        minimum_singular_value = float(min(cross_geometry_singular_values))
        maximum_singular_value = float(max(cross_geometry_singular_values))
    else:
        minimum_singular_value = None
        maximum_singular_value = None
    return MolecularSOCReplayOverlapDiagnosticsV232(
        record_count=int(nrecord_left),
        state_count=int(nstate_left),
        unordered_cross_geometry_pair_count=int(
            nrecord_left * (nrecord_left - 1) // 2
        ),
        maximum_self_identity_residual=float(max(self_residuals, default=0.0)),
        maximum_reciprocity_residual=float(
            max(reciprocity_residuals, default=0.0)
        ),
        minimum_cross_geometry_singular_value=minimum_singular_value,
        maximum_cross_geometry_singular_value=maximum_singular_value,
        maximum_contraction_excess=float(
            max(0.0, (maximum_singular_value or 0.0) - 1.0)
        ),
    )


@dataclass(frozen=True)
class MolecularSOCReplayDatasetV230:
    q: np.ndarray
    H_spin_free: np.ndarray
    K_spin_free: np.ndarray
    H_soc: np.ndarray
    K_soc: np.ndarray
    connection_q: np.ndarray
    mass_matrix_q_au: np.ndarray
    overlaps: np.ndarray
    converged: np.ndarray
    provenance: ElectronicOperatorProvenanceV213
    symmetry_contract: SOCSymmetryContractV221
    molecular_soc_contract: MolecularSOCAdmissionContractV230
    dataset_fingerprint: str
    manifest_path: Path
    arrays_path: Path
    coordinate_digits: int = 14

    def overlap_diagnostics(self):
        """Return v0.23.2 diagnostics without assuming cross-geometry unitarity."""
        return _replay_overlap_diagnostics_v232(self.overlaps)

    def validate(self, *, tolerance=1.0e-10):
        tolerance = float(tolerance)
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("replay validation tolerance must be finite and nonnegative.")
        q = np.asarray(self.q, dtype=float)
        H0 = np.asarray(self.H_spin_free, dtype=complex)
        K0 = np.asarray(self.K_spin_free, dtype=complex)
        Hso = np.asarray(self.H_soc, dtype=complex)
        Kso = np.asarray(self.K_soc, dtype=complex)
        D = np.asarray(self.connection_q, dtype=complex)
        mass = np.asarray(self.mass_matrix_q_au, dtype=float)
        overlaps = np.asarray(self.overlaps, dtype=complex)
        converged = np.asarray(self.converged)
        if q.ndim != 2 or q.shape[0] < 1 or q.shape[1] < 1:
            raise ValueError("replay coordinates must have shape (nrecord,nq).")
        nrecord, nq = q.shape
        if H0.ndim != 3 or H0.shape[0] != nrecord or H0.shape[1] != H0.shape[2]:
            raise ValueError("replay spin-free Hamiltonians have incompatible shape.")
        nstate = H0.shape[1]
        if nstate < 1:
            raise ValueError("replay state dimension must be positive.")
        expected_H = (nrecord, nstate, nstate)
        expected_K = (nrecord, nq, nstate, nstate)
        if Hso.shape != expected_H:
            raise ValueError("replay SOC Hamiltonians have incompatible shape.")
        if K0.shape != expected_K or Kso.shape != expected_K or D.shape != expected_K:
            raise ValueError("replay derivative arrays have incompatible shape.")
        if mass.shape != (nrecord, nq, nq):
            raise ValueError("replay mass matrices have incompatible shape.")
        if overlaps.shape != (nrecord, nrecord, nstate, nstate):
            raise ValueError("replay overlap table has incompatible shape.")
        if converged.shape != (nrecord,) or converged.dtype.kind != "b":
            raise ValueError("replay convergence flags must be a Boolean record vector.")
        for array, name in (
            (q, "q"),
            (H0, "H_spin_free"),
            (K0, "K_spin_free"),
            (Hso, "H_soc"),
            (Kso, "K_soc"),
            (D, "connection_q"),
            (mass, "mass_matrix_q_au"),
            (overlaps, "overlaps"),
        ):
            if not np.all(np.isfinite(array)):
                raise ValueError(f"replay {name} contains non-finite data.")
        keys = [tuple(np.round(row, int(self.coordinate_digits))) for row in q]
        if len(set(keys)) != nrecord:
            raise ValueError("replay coordinates collide under the coordinate policy.")
        overlap_diagnostics = _replay_overlap_diagnostics_v232(overlaps)
        if overlap_diagnostics.maximum_self_identity_residual > tolerance:
            raise ValueError(
                "replay self overlaps violate the exact identity contract."
            )
        if overlap_diagnostics.maximum_reciprocity_residual > tolerance:
            raise ValueError(
                "replay cross-geometry overlaps violate adjoint reciprocity."
            )
        if overlap_diagnostics.maximum_contraction_excess > tolerance:
            raise ValueError(
                "replay cross-geometry overlap is expansive; singular values "
                "must not exceed one."
            )
        if self.provenance.model_space.nstate != nstate:
            raise ValueError("replay state dimension disagrees with provenance.")
        self.molecular_soc_contract.validate(self.symmetry_contract)
        if self.molecular_soc_contract.all_electronic_calculations_converged != bool(
            np.all(converged)
        ):
            raise ValueError(
                "replay record convergence flags disagree with the admission contract."
            )
        parameters = self.provenance.parameters
        if parameters.get("v230_molecular_soc_contract_fingerprint") != (
            self.molecular_soc_contract.fingerprint()
        ):
            raise ValueError("replay provenance and molecular SOC contract disagree.")
        if not str(self.dataset_fingerprint).strip():
            raise ValueError("replay dataset fingerprint cannot be empty.")
        return self


def capture_molecular_soc_replay_v230(
    directory,
    source_provider,
    coordinates,
    molecular_soc_contract,
    *,
    provenance=None,
    convergence_flags=None,
    coordinate_digits=14,
    overwrite=False,
):
    """Capture exact provider records and all pair overlaps into a replay dataset."""
    directory = Path(directory)
    coordinate_digits = int(coordinate_digits)
    if coordinate_digits < 0:
        raise ValueError("coordinate_digits cannot be negative.")
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[0] < 1 or coordinates.shape[1] < 1:
        raise ValueError("coordinates must have shape (nrecord,nq).")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("replay coordinates must be finite.")
    symmetry = getattr(source_provider, "soc_symmetry_contract", None)
    if callable(symmetry):
        symmetry = symmetry()
    if not isinstance(symmetry, SOCSymmetryContractV221):
        raise TypeError("source provider lacks a v0.22.1 SOC symmetry contract.")
    molecular_soc_contract = molecular_soc_contract.validate(symmetry)
    if not molecular_soc_contract.capabilities.deterministic_replay:
        raise ValueError("replay capture requires deterministic_replay capability.")
    if provenance is None:
        provenance = provenance_with_molecular_soc_contract_v230(
            source_provider.provenance, molecular_soc_contract
        )
    provenance = provenance.validate()

    snapshots = []
    components = []
    for q in coordinates:
        components.append(source_provider.components(q).validate())
        snapshots.append(source_provider.evaluate_snapshot(q).validate())
    nrecord = len(snapshots)
    nstate = components[0].H.shape[0]
    nq = coordinates.shape[1]
    if convergence_flags is None:
        if molecular_soc_contract.identity.source_kind != "validation_fixture":
            raise ValueError(
                "real molecular replay capture requires explicit per-record convergence flags."
            )
        convergence_flags = np.ones(nrecord, dtype=bool)
    convergence_flags = np.asarray(convergence_flags)
    if convergence_flags.shape != (nrecord,) or convergence_flags.dtype.kind != "b":
        raise ValueError("convergence_flags must be a Boolean value per record.")
    if bool(np.all(convergence_flags)) != bool(
        molecular_soc_contract.all_electronic_calculations_converged
    ):
        raise ValueError("convergence_flags disagree with the admission contract.")

    overlaps = np.empty((nrecord, nrecord, nstate, nstate), dtype=complex)
    for left in range(nrecord):
        for right in range(nrecord):
            overlaps[left, right] = np.asarray(
                source_provider.snapshot_overlap(snapshots[left], snapshots[right]),
                dtype=complex,
            )
    projector_names = tuple(sorted(symmetry.projectors))
    arrays = {
        "q": coordinates,
        "H_spin_free": np.stack([item.H_spin_free for item in components]),
        "K_spin_free": np.stack([item.K_spin_free for item in components]),
        "H_soc": np.stack([item.H_soc for item in components]),
        "K_soc": np.stack([item.K_soc for item in components]),
        "connection_q": np.stack(
            [snapshot.point.connection_q for snapshot in snapshots]
        ),
        "mass_matrix_q_au": np.stack(
            [snapshot.point.mass_matrix_q_au for snapshot in snapshots]
        ),
        "overlaps": overlaps,
        "converged": convergence_flags.astype(bool),
        "time_reversal_matrix": symmetry.time_reversal_matrix,
        "projectors": np.stack([symmetry.projectors[name] for name in projector_names]),
    }
    if any(array.shape[0] != nrecord for name, array in arrays.items() if name not in {"time_reversal_matrix", "projectors"}):
        raise ValueError("source provider changed dimension during replay capture.")

    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory / REPLAY_MANIFEST_NAME_V230
    arrays_path = directory / REPLAY_ARRAYS_NAME_V230
    if not overwrite and (manifest_path.exists() or arrays_path.exists()):
        raise FileExistsError("molecular SOC replay target already exists.")
    temporary_arrays = directory / f".{REPLAY_ARRAYS_NAME_V230}.tmp"
    temporary_manifest = directory / f".{REPLAY_MANIFEST_NAME_V230}.tmp"
    _write_deterministic_npz_v230(temporary_arrays, arrays)
    arrays_sha256 = _sha256_file_v230(temporary_arrays)
    payload = {
        "format": "gaussian-nadyn-molecular-soc-replay",
        "format_version": 1,
        "release": "v0.23.0",
        "arrays_file": REPLAY_ARRAYS_NAME_V230,
        "arrays_sha256": arrays_sha256,
        "coordinate_digits": coordinate_digits,
        "nrecord": nrecord,
        "nq": nq,
        "nstate": nstate,
        "projector_names": list(projector_names),
        "electron_parity": symmetry.electron_parity,
        "external_magnetic_field": bool(symmetry.external_magnetic_field),
        "provenance": provenance.as_dict(),
        "provenance_fingerprint": provenance.fingerprint(),
        "molecular_soc_contract": molecular_soc_contract.as_dict(),
        "molecular_soc_contract_fingerprint": molecular_soc_contract.fingerprint(),
        "source_provider_fingerprint": source_provider.provenance.fingerprint(),
    }
    dataset_fingerprint = hashlib.sha256(_canonical_json_bytes_v230(payload)).hexdigest()
    manifest = {**payload, "dataset_fingerprint": dataset_fingerprint}
    temporary_manifest.write_bytes(_canonical_json_bytes_v230(manifest))
    os.replace(temporary_arrays, arrays_path)
    os.replace(temporary_manifest, manifest_path)
    return load_molecular_soc_replay_v230(manifest_path)


def load_molecular_soc_replay_v230(manifest_path):
    manifest_path = Path(manifest_path)
    if manifest_path.is_dir():
        manifest_path = manifest_path / REPLAY_MANIFEST_NAME_V230
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != "gaussian-nadyn-molecular-soc-replay":
        raise ValueError("molecular SOC replay format mismatch.")
    if manifest.get("format_version") != 1:
        raise ValueError("molecular SOC replay version mismatch.")
    arrays_path = manifest_path.parent / manifest["arrays_file"]
    if _sha256_file_v230(arrays_path) != manifest["arrays_sha256"]:
        raise ValueError("molecular SOC replay array integrity check failed.")
    fingerprint_payload = dict(manifest)
    stored_dataset_fingerprint = fingerprint_payload.pop("dataset_fingerprint", None)
    computed_dataset_fingerprint = hashlib.sha256(
        _canonical_json_bytes_v230(fingerprint_payload)
    ).hexdigest()
    if stored_dataset_fingerprint != computed_dataset_fingerprint:
        raise ValueError("molecular SOC replay manifest integrity check failed.")
    provenance = _provenance_from_dict_v230(manifest["provenance"])
    if provenance.fingerprint() != manifest["provenance_fingerprint"]:
        raise ValueError("molecular SOC replay provenance fingerprint mismatch.")
    contract = _contract_from_dict_v230(manifest["molecular_soc_contract"])
    if contract.fingerprint() != manifest["molecular_soc_contract_fingerprint"]:
        raise ValueError("molecular SOC replay contract fingerprint mismatch.")
    with np.load(arrays_path, allow_pickle=False) as arrays:
        required = {
            "q",
            "H_spin_free",
            "K_spin_free",
            "H_soc",
            "K_soc",
            "connection_q",
            "mass_matrix_q_au",
            "overlaps",
            "converged",
            "time_reversal_matrix",
            "projectors",
        }
        if set(arrays.files) != required:
            raise ValueError("molecular SOC replay array member set mismatch.")
        data = {name: arrays[name].copy() for name in required}
    names = tuple(manifest["projector_names"])
    if data["projectors"].shape[0] != len(names) or len(set(names)) != len(names):
        raise ValueError("molecular SOC replay projector names are inconsistent.")
    symmetry = SOCSymmetryContractV221(
        electron_parity=manifest["electron_parity"],
        time_reversal_matrix=data["time_reversal_matrix"],
        projectors={
            name: data["projectors"][index] for index, name in enumerate(names)
        },
        external_magnetic_field=manifest["external_magnetic_field"],
    )
    if manifest["nrecord"] != data["q"].shape[0]:
        raise ValueError("molecular SOC replay record count mismatch.")
    if manifest["nq"] != data["q"].shape[1]:
        raise ValueError("molecular SOC replay coordinate dimension mismatch.")
    if manifest["nstate"] != data["H_spin_free"].shape[1]:
        raise ValueError("molecular SOC replay state dimension mismatch.")
    return MolecularSOCReplayDatasetV230(
        q=data["q"],
        H_spin_free=data["H_spin_free"],
        K_spin_free=data["K_spin_free"],
        H_soc=data["H_soc"],
        K_soc=data["K_soc"],
        connection_q=data["connection_q"],
        mass_matrix_q_au=data["mass_matrix_q_au"],
        overlaps=data["overlaps"],
        converged=data["converged"],
        provenance=provenance,
        symmetry_contract=symmetry,
        molecular_soc_contract=contract,
        dataset_fingerprint=stored_dataset_fingerprint,
        manifest_path=manifest_path.resolve(),
        arrays_path=arrays_path.resolve(),
        coordinate_digits=int(manifest["coordinate_digits"]),
    ).validate()


class FileBackedMolecularSOCProviderV230:
    """Exact-record replay provider; interpolation and extrapolation are forbidden."""

    def __init__(self, manifest_path):
        self.dataset = load_molecular_soc_replay_v230(manifest_path)
        self.provenance = self.dataset.provenance
        self._keys = {
            tuple(np.round(row, self.dataset.coordinate_digits)): index
            for index, row in enumerate(self.dataset.q)
        }

    @property
    def molecular_soc_contract(self):
        return self.dataset.molecular_soc_contract

    @property
    def soc_symmetry_contract(self):
        return self.dataset.symmetry_contract

    @property
    def time_reversal_matrix(self):
        return self.dataset.symmetry_contract.time_reversal_matrix.copy()

    @property
    def projectors(self):
        return {
            name: value.copy()
            for name, value in self.dataset.symmetry_contract.projectors.items()
        }

    @property
    def replay_fingerprint(self):
        return self.dataset.dataset_fingerprint

    def _index(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (self.dataset.q.shape[1],) or not np.all(np.isfinite(q)):
            raise ValueError("molecular SOC replay request has incompatible coordinates.")
        key = tuple(np.round(q, self.dataset.coordinate_digits))
        index = self._keys.get(key)
        if index is None or not np.allclose(
            q, self.dataset.q[index], rtol=0.0, atol=10.0 ** (-self.dataset.coordinate_digits)
        ):
            raise KeyError(
                "molecular SOC replay contains no exact record for the requested geometry; "
                "interpolation and extrapolation are forbidden."
            )
        return index

    def components(self, q):
        index = self._index(q)
        return SOCOperatorComponentsV220(
            self.dataset.q[index].copy(),
            self.dataset.H_spin_free[index].copy(),
            self.dataset.K_spin_free[index].copy(),
            self.dataset.H_soc[index].copy(),
            self.dataset.K_soc[index].copy(),
        ).validate()

    def evaluate_snapshot(self, q):
        index = self._index(q)
        components = self.components(q)
        point = compose_electronic_operator_v213(
            q=components.q,
            H_spin_free=components.H_spin_free,
            dH_spin_free_dq=components.K_spin_free,
            H_soc=components.H_soc,
            dH_soc_dq=components.K_soc,
            connection_q=self.dataset.connection_q[index],
            mass_matrix_q_au=self.dataset.mass_matrix_q_au[index],
            provenance=self.provenance,
        )
        point.metadata.update(
            {
                "v230_replay_dataset_fingerprint": self.replay_fingerprint,
                "v230_replay_record_index": int(index),
                "v230_electronic_converged": bool(self.dataset.converged[index]),
            }
        )
        return ElectronicOperatorSnapshotV21(
            point=point,
            wavefunction_snapshot=ReplayWavefunctionTokenV230(
                self.replay_fingerprint, int(index)
            ),
            metadata={
                "provider": "FileBackedMolecularSOCProviderV230",
                "dataset_fingerprint": self.replay_fingerprint,
                "record_index": int(index),
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def snapshot_overlap(self, left, right):
        left_token = left.wavefunction_snapshot
        right_token = right.wavefunction_snapshot
        if not isinstance(left_token, ReplayWavefunctionTokenV230) or not isinstance(
            right_token, ReplayWavefunctionTokenV230
        ):
            raise TypeError("replay overlap requires replay snapshot tokens.")
        if (
            left_token.dataset_fingerprint != self.replay_fingerprint
            or right_token.dataset_fingerprint != self.replay_fingerprint
        ):
            raise ValueError("cross-dataset molecular SOC overlap is forbidden.")
        return self.dataset.overlaps[
            left_token.record_index, right_token.record_index
        ].copy()

    def diagnostics_dict(self):
        overlap_diagnostics = self.dataset.overlap_diagnostics()
        return {
            "provider": "FileBackedMolecularSOCProviderV230",
            "dataset_fingerprint": self.replay_fingerprint,
            "records": int(len(self.dataset.q)),
            "exact_record_only": True,
            "all_electronic_calculations_converged": bool(
                np.all(self.dataset.converged)
            ),
            "overlap_contract_v232": overlap_diagnostics.as_dict(),
        }
