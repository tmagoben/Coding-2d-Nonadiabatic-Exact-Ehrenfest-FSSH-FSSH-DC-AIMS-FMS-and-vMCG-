"""Versioned v0.23.3 replay with certified overlap transports and NAC identity."""

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .analytic_soc_models_v220 import SOCOperatorComponentsV220
from .electronic_contract_v213 import compose_electronic_operator_v213
from .electronic_operator_v21 import ElectronicOperatorSnapshotV21
from .finite_manifold_transport_v233 import (
    OVERLAP_CONTRACT_ID_V233,
    TRANSPORT_CONTRACT_ID_V233,
    FiniteManifoldOverlapPolicyV233,
    analyze_finite_manifold_overlap_v233,
    certify_reciprocal_transport_pair_v233,
)
from .molecular_soc_replay_v230 import (
    MolecularSOCReplayDatasetV230,
    _canonical_json_bytes_v230,
    _contract_from_dict_v230,
    _provenance_from_dict_v230,
    _sha256_file_v230,
    _write_deterministic_npz_v230,
    capture_molecular_soc_replay_v230,
    load_molecular_soc_replay_v230,
)
from .nac_compatibility_v233 import (
    DerivativeCouplingConventionV233,
    LegacyReplayMigrationAttestationV233,
    derivative_coupling_convention_from_dict_v233,
)
from .provider_numerical_identity_v233 import (
    build_provider_numerical_identity_v233,
)
from .soc_admission_v221 import SOCSymmetryContractV221


REPLAY_MANIFEST_NAME_V233 = "molecular_soc_manifest_v233.json"
REPLAY_ARRAYS_NAME_V233 = "molecular_soc_arrays_v233.npz"
REPLAY_FORMAT_V233 = "gaussian-nadyn-molecular-soc-replay"
REPLAY_FORMAT_VERSION_V233 = 2


def _policy_from_dict_v233(payload):
    if not isinstance(payload, dict):
        raise TypeError("overlap-quality policy payload must be a mapping.")
    expected = set(FiniteManifoldOverlapPolicyV233().__dict__)
    if set(payload) != expected:
        raise ValueError("overlap-quality policy field set mismatch.")
    return FiniteManifoldOverlapPolicyV233(**payload).validate()


def _symmetry_from_arrays_v233(manifest, arrays):
    names = tuple(manifest["projector_names"])
    projectors = arrays["projectors"]
    if projectors.shape[0] != len(names) or len(set(names)) != len(names):
        raise ValueError("v0.23.3 replay projector names are inconsistent.")
    return SOCSymmetryContractV221(
        electron_parity=manifest["electron_parity"],
        time_reversal_matrix=arrays["time_reversal_matrix"],
        projectors={name: projectors[index] for index, name in enumerate(names)},
        external_magnetic_field=manifest["external_magnetic_field"],
    )


@dataclass(frozen=True)
class ReplayWavefunctionTokenV233:
    dataset_fingerprint: str
    record_index: int


@dataclass(frozen=True)
class MolecularSOCReplayDatasetV233:
    legacy: MolecularSOCReplayDatasetV230
    overlap_transports: np.ndarray
    overlap_singular_values: np.ndarray
    overlap_policy: FiniteManifoldOverlapPolicyV233
    nac_convention: DerivativeCouplingConventionV233
    dataset_fingerprint: str
    legacy_dataset_fingerprint: str
    manifest_path: Path
    arrays_path: Path
    migration_attestation: dict | None = None

    @property
    def q(self):
        return self.legacy.q

    @property
    def H_spin_free(self):
        return self.legacy.H_spin_free

    @property
    def K_spin_free(self):
        return self.legacy.K_spin_free

    @property
    def H_soc(self):
        return self.legacy.H_soc

    @property
    def K_soc(self):
        return self.legacy.K_soc

    @property
    def connection_q(self):
        return self.legacy.connection_q

    @property
    def mass_matrix_q_au(self):
        return self.legacy.mass_matrix_q_au

    @property
    def overlaps(self):
        return self.legacy.overlaps

    @property
    def converged(self):
        return self.legacy.converged

    @property
    def provenance(self):
        return self.legacy.provenance

    @property
    def symmetry_contract(self):
        return self.legacy.symmetry_contract

    @property
    def molecular_soc_contract(self):
        return self.legacy.molecular_soc_contract

    @property
    def coordinate_digits(self):
        return self.legacy.coordinate_digits

    def overlap_diagnostics(self):
        return self.legacy.overlap_diagnostics()

    def transport(self, left_index, right_index):
        return np.asarray(
            self.overlap_transports[int(left_index), int(right_index)],
            dtype=complex,
        ).copy()

    def validate(self, *, tolerance=1.0e-10):
        self.legacy.validate(tolerance=tolerance)
        policy = self.overlap_policy.validate()
        convention = self.nac_convention.validate()
        if convention.use_etfs:
            raise ValueError(
                "trajectory-ready replay requires full-overlap NACs with "
                "use_etfs=False."
            )
        nrecord = len(self.q)
        nstate = self.H_spin_free.shape[1]
        transports = np.asarray(self.overlap_transports, dtype=complex)
        singular_values = np.asarray(self.overlap_singular_values, dtype=float)
        if transports.shape != (nrecord, nrecord, nstate, nstate):
            raise ValueError("v0.23.3 replay transport table has incompatible shape.")
        if singular_values.shape != (nrecord, nrecord, nstate):
            raise ValueError(
                "v0.23.3 replay singular-value table has incompatible shape."
            )
        if not np.all(np.isfinite(transports)) or not np.all(
            np.isfinite(singular_values)
        ):
            raise ValueError("v0.23.3 replay transport data are non-finite.")
        identity = np.eye(nstate, dtype=complex)
        for index in range(nrecord):
            if np.linalg.norm(transports[index, index] - identity, ord="fro") > tolerance:
                raise ValueError("v0.23.3 self transport is not identity.")
            if np.linalg.norm(singular_values[index, index] - 1.0) > tolerance:
                raise ValueError("v0.23.3 self-overlap singular values are not one.")
        for left in range(nrecord):
            for right in range(left + 1, nrecord):
                pair = certify_reciprocal_transport_pair_v233(
                    self.overlaps[left, right],
                    self.overlaps[right, left],
                    policy=policy,
                    reciprocity_tolerance=tolerance,
                )
                for row, col, result in (
                    (left, right, pair.left_to_right_block),
                    (right, left, pair.right_to_left_block),
                ):
                    if np.linalg.norm(
                        transports[row, col]
                        - result.right_to_left_transport,
                        ord="fro",
                    ) > tolerance:
                        raise ValueError(
                            "stored v0.23.3 transport differs from the raw-overlap "
                            "polar factor."
                        )
                    if np.linalg.norm(
                        singular_values[row, col] - result.singular_values
                    ) > tolerance:
                        raise ValueError(
                            "stored v0.23.3 singular values differ from the raw overlap."
                        )
        if not str(self.dataset_fingerprint).strip():
            raise ValueError("v0.23.3 replay fingerprint cannot be empty.")
        if not str(self.legacy_dataset_fingerprint).strip():
            raise ValueError("v0.23.3 legacy-source fingerprint cannot be empty.")
        return self


def _write_v233_from_legacy_dataset(
    directory,
    legacy,
    *,
    nac_convention,
    overlap_policy,
    migration_attestation=None,
    overwrite=False,
):
    directory = Path(directory)
    legacy = legacy.validate()
    if type(nac_convention) is not DerivativeCouplingConventionV233:
        raise TypeError("nac_convention must be DerivativeCouplingConventionV233.")
    nac_convention = nac_convention.validate()
    if nac_convention.use_etfs:
        raise ValueError("v0.23.3 replay capture requires use_etfs=False.")
    overlap_policy = overlap_policy.validate()
    nrecord = len(legacy.q)
    nstate = legacy.H_spin_free.shape[1]
    transports = np.empty((nrecord, nrecord, nstate, nstate), dtype=complex)
    singular_values = np.empty((nrecord, nrecord, nstate), dtype=float)
    identity = np.eye(nstate, dtype=complex)
    for index in range(nrecord):
        transports[index, index] = identity
        singular_values[index, index] = 1.0
    for left in range(nrecord):
        for right in range(left + 1, nrecord):
            pair = certify_reciprocal_transport_pair_v233(
                legacy.overlaps[left, right],
                legacy.overlaps[right, left],
                policy=overlap_policy,
            )
            transports[left, right] = (
                pair.left_to_right_block.right_to_left_transport
            )
            transports[right, left] = (
                pair.right_to_left_block.right_to_left_transport
            )
            singular_values[left, right] = pair.left_to_right_block.singular_values
            singular_values[right, left] = pair.right_to_left_block.singular_values

    symmetry = legacy.symmetry_contract
    projector_names = tuple(sorted(symmetry.projectors))
    arrays = {
        "q": legacy.q,
        "H_spin_free": legacy.H_spin_free,
        "K_spin_free": legacy.K_spin_free,
        "H_soc": legacy.H_soc,
        "K_soc": legacy.K_soc,
        "connection_q": legacy.connection_q,
        "mass_matrix_q_au": legacy.mass_matrix_q_au,
        "overlaps": legacy.overlaps,
        "overlap_transports": transports,
        "overlap_singular_values": singular_values,
        "converged": legacy.converged,
        "time_reversal_matrix": symmetry.time_reversal_matrix,
        "projectors": np.stack(
            [symmetry.projectors[name] for name in projector_names]
        ),
    }
    directory.mkdir(parents=True, exist_ok=True)
    manifest_path = directory / REPLAY_MANIFEST_NAME_V233
    arrays_path = directory / REPLAY_ARRAYS_NAME_V233
    if not overwrite and (manifest_path.exists() or arrays_path.exists()):
        raise FileExistsError("v0.23.3 molecular SOC replay target already exists.")
    temporary_arrays = directory / f".{REPLAY_ARRAYS_NAME_V233}.tmp"
    temporary_manifest = directory / f".{REPLAY_MANIFEST_NAME_V233}.tmp"
    _write_deterministic_npz_v230(temporary_arrays, arrays)
    arrays_sha256 = _sha256_file_v230(temporary_arrays)
    payload = {
        "format": REPLAY_FORMAT_V233,
        "format_version": REPLAY_FORMAT_VERSION_V233,
        "release": "v0.23.3",
        "arrays_file": REPLAY_ARRAYS_NAME_V233,
        "arrays_sha256": arrays_sha256,
        "coordinate_digits": int(legacy.coordinate_digits),
        "nrecord": int(nrecord),
        "nq": int(legacy.q.shape[1]),
        "nstate": int(nstate),
        "projector_names": list(projector_names),
        "electron_parity": symmetry.electron_parity,
        "external_magnetic_field": bool(symmetry.external_magnetic_field),
        "provenance": legacy.provenance.as_dict(),
        "provenance_fingerprint": legacy.provenance.fingerprint(),
        "molecular_soc_contract": legacy.molecular_soc_contract.as_dict(),
        "molecular_soc_contract_fingerprint": (
            legacy.molecular_soc_contract.fingerprint()
        ),
        "legacy_dataset_fingerprint": legacy.dataset_fingerprint,
        "overlap_contract": OVERLAP_CONTRACT_ID_V233,
        "transport_contract": TRANSPORT_CONTRACT_ID_V233,
        "overlap_quality_policy": overlap_policy.as_dict(),
        "nac_convention": nac_convention.as_dict(),
        "nac_convention_fingerprint": nac_convention.fingerprint(),
        "migration_attestation": migration_attestation,
    }
    dataset_fingerprint = hashlib.sha256(
        _canonical_json_bytes_v230(payload)
    ).hexdigest()
    manifest = {**payload, "dataset_fingerprint": dataset_fingerprint}
    temporary_manifest.write_bytes(_canonical_json_bytes_v230(manifest))
    os.replace(temporary_arrays, arrays_path)
    os.replace(temporary_manifest, manifest_path)
    return load_molecular_soc_replay_v233(manifest_path)


def capture_molecular_soc_replay_v233(
    directory,
    source_provider,
    coordinates,
    molecular_soc_contract,
    *,
    nac_convention,
    overlap_policy=FiniteManifoldOverlapPolicyV233(),
    provenance=None,
    convergence_flags=None,
    coordinate_digits=14,
    overwrite=False,
):
    """Capture new data directly into the v0.23.3 replay semantics."""
    with tempfile.TemporaryDirectory(prefix="gnd-v233-capture-") as temporary:
        legacy = capture_molecular_soc_replay_v230(
            temporary,
            source_provider,
            coordinates,
            molecular_soc_contract,
            provenance=provenance,
            convergence_flags=convergence_flags,
            coordinate_digits=coordinate_digits,
        )
        return _write_v233_from_legacy_dataset(
            directory,
            legacy,
            nac_convention=nac_convention,
            overlap_policy=overlap_policy,
            migration_attestation=None,
            overwrite=overwrite,
        )


def migrate_molecular_soc_replay_v230_to_v233(
    legacy_manifest_path,
    output_directory,
    *,
    nac_convention,
    migration_attestation,
    overlap_policy=FiniteManifoldOverlapPolicyV233(),
    overwrite=False,
):
    """Migrate only explicitly attested legacy data; never infer an NAC sign."""
    legacy = load_molecular_soc_replay_v230(legacy_manifest_path)
    if type(migration_attestation) is not LegacyReplayMigrationAttestationV233:
        raise TypeError(
            "migration_attestation must be LegacyReplayMigrationAttestationV233."
        )
    migration_attestation.validate(
        expected_legacy_fingerprint=legacy.dataset_fingerprint,
        convention=nac_convention,
    )
    return _write_v233_from_legacy_dataset(
        output_directory,
        legacy,
        nac_convention=nac_convention,
        overlap_policy=overlap_policy,
        migration_attestation=migration_attestation.as_dict(),
        overwrite=overwrite,
    )


def load_molecular_soc_replay_v233(manifest_path):
    manifest_path = Path(manifest_path)
    if manifest_path.is_dir():
        manifest_path = manifest_path / REPLAY_MANIFEST_NAME_V233
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != REPLAY_FORMAT_V233:
        raise ValueError("v0.23.3 molecular SOC replay format mismatch.")
    if manifest.get("format_version") != REPLAY_FORMAT_VERSION_V233:
        if manifest.get("format_version") == 1:
            raise ValueError(
                "legacy v0.23.0 replay is quarantined; use the explicit v0.23.3 "
                "migration API with an NAC attestation."
            )
        raise ValueError("v0.23.3 molecular SOC replay version mismatch.")
    if manifest.get("overlap_contract") != OVERLAP_CONTRACT_ID_V233:
        raise ValueError("v0.23.3 replay overlap contract mismatch.")
    if manifest.get("transport_contract") != TRANSPORT_CONTRACT_ID_V233:
        raise ValueError("v0.23.3 replay transport contract mismatch.")
    arrays_path = manifest_path.parent / manifest["arrays_file"]
    if _sha256_file_v230(arrays_path) != manifest["arrays_sha256"]:
        raise ValueError("v0.23.3 replay array integrity check failed.")
    fingerprint_payload = dict(manifest)
    stored_dataset_fingerprint = fingerprint_payload.pop(
        "dataset_fingerprint", None
    )
    computed_dataset_fingerprint = hashlib.sha256(
        _canonical_json_bytes_v230(fingerprint_payload)
    ).hexdigest()
    if stored_dataset_fingerprint != computed_dataset_fingerprint:
        raise ValueError("v0.23.3 replay manifest integrity check failed.")
    provenance = _provenance_from_dict_v230(manifest["provenance"])
    if provenance.fingerprint() != manifest["provenance_fingerprint"]:
        raise ValueError("v0.23.3 replay provenance fingerprint mismatch.")
    contract = _contract_from_dict_v230(manifest["molecular_soc_contract"])
    if contract.fingerprint() != manifest["molecular_soc_contract_fingerprint"]:
        raise ValueError("v0.23.3 replay contract fingerprint mismatch.")
    convention = derivative_coupling_convention_from_dict_v233(
        manifest["nac_convention"]
    )
    if convention.fingerprint() != manifest["nac_convention_fingerprint"]:
        raise ValueError("v0.23.3 replay NAC convention fingerprint mismatch.")
    policy = _policy_from_dict_v233(manifest["overlap_quality_policy"])
    required = {
        "q",
        "H_spin_free",
        "K_spin_free",
        "H_soc",
        "K_soc",
        "connection_q",
        "mass_matrix_q_au",
        "overlaps",
        "overlap_transports",
        "overlap_singular_values",
        "converged",
        "time_reversal_matrix",
        "projectors",
    }
    with np.load(arrays_path, allow_pickle=False) as archive:
        if set(archive.files) != required:
            raise ValueError("v0.23.3 replay array member set mismatch.")
        arrays = {name: archive[name].copy() for name in required}
    symmetry = _symmetry_from_arrays_v233(manifest, arrays)
    if manifest["nrecord"] != arrays["q"].shape[0]:
        raise ValueError("v0.23.3 replay record count mismatch.")
    if manifest["nq"] != arrays["q"].shape[1]:
        raise ValueError("v0.23.3 replay coordinate dimension mismatch.")
    if manifest["nstate"] != arrays["H_spin_free"].shape[1]:
        raise ValueError("v0.23.3 replay state dimension mismatch.")
    embedded_legacy = MolecularSOCReplayDatasetV230(
        q=arrays["q"],
        H_spin_free=arrays["H_spin_free"],
        K_spin_free=arrays["K_spin_free"],
        H_soc=arrays["H_soc"],
        K_soc=arrays["K_soc"],
        connection_q=arrays["connection_q"],
        mass_matrix_q_au=arrays["mass_matrix_q_au"],
        overlaps=arrays["overlaps"],
        converged=arrays["converged"],
        provenance=provenance,
        symmetry_contract=symmetry,
        molecular_soc_contract=contract,
        dataset_fingerprint=manifest["legacy_dataset_fingerprint"],
        manifest_path=manifest_path.resolve(),
        arrays_path=arrays_path.resolve(),
        coordinate_digits=int(manifest["coordinate_digits"]),
    )
    return MolecularSOCReplayDatasetV233(
        legacy=embedded_legacy,
        overlap_transports=arrays["overlap_transports"],
        overlap_singular_values=arrays["overlap_singular_values"],
        overlap_policy=policy,
        nac_convention=convention,
        dataset_fingerprint=stored_dataset_fingerprint,
        legacy_dataset_fingerprint=manifest["legacy_dataset_fingerprint"],
        manifest_path=manifest_path.resolve(),
        arrays_path=arrays_path.resolve(),
        migration_attestation=manifest.get("migration_attestation"),
    ).validate()


class FileBackedMolecularSOCProviderV233:
    """Exact v0.23.3 replay exposing raw overlap and certified transport separately."""

    def __init__(self, manifest_path):
        self.dataset = load_molecular_soc_replay_v233(manifest_path)
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
    def replay_fingerprint(self):
        return self.dataset.dataset_fingerprint

    @property
    def numerical_identity_v233(self):
        return build_provider_numerical_identity_v233(
            self.provenance,
            self.dataset.nac_convention,
            overlap_policy=self.dataset.overlap_policy,
        )

    @property
    def time_reversal_matrix(self):
        return self.dataset.symmetry_contract.time_reversal_matrix.copy()

    @property
    def projectors(self):
        return {
            name: value.copy()
            for name, value in self.dataset.symmetry_contract.projectors.items()
        }

    def _index(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (self.dataset.q.shape[1],) or not np.all(np.isfinite(q)):
            raise ValueError("v0.23.3 replay request has incompatible coordinates.")
        key = tuple(np.round(q, self.dataset.coordinate_digits))
        index = self._keys.get(key)
        if index is None or not np.allclose(
            q,
            self.dataset.q[index],
            rtol=0.0,
            atol=10.0 ** (-self.dataset.coordinate_digits),
        ):
            raise KeyError(
                "v0.23.3 replay contains no exact record for the requested geometry."
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
                "v233_replay_dataset_fingerprint": self.replay_fingerprint,
                "v233_replay_record_index": int(index),
                "v233_electronic_converged": bool(
                    self.dataset.converged[index]
                ),
                "v233_nac_convention_fingerprint": (
                    self.dataset.nac_convention.fingerprint()
                ),
                "v233_provider_numerical_identity_fingerprint": (
                    self.numerical_identity_v233.fingerprint()
                ),
            }
        )
        return ElectronicOperatorSnapshotV21(
            point=point,
            wavefunction_snapshot=ReplayWavefunctionTokenV233(
                self.replay_fingerprint, int(index)
            ),
            metadata={
                "provider": "FileBackedMolecularSOCProviderV233",
                "dataset_fingerprint": self.replay_fingerprint,
                "record_index": int(index),
                "v233_nac_convention_fingerprint": (
                    self.dataset.nac_convention.fingerprint()
                ),
                "v233_provider_numerical_identity_fingerprint": (
                    self.numerical_identity_v233.fingerprint()
                ),
            },
        ).validate()

    def evaluate(self, q):
        return self.evaluate_snapshot(q).point

    def _tokens(self, left, right):
        left_token = left.wavefunction_snapshot
        right_token = right.wavefunction_snapshot
        if not isinstance(left_token, ReplayWavefunctionTokenV233) or not isinstance(
            right_token, ReplayWavefunctionTokenV233
        ):
            raise TypeError("v0.23.3 replay requires v0.23.3 snapshot tokens.")
        if (
            left_token.dataset_fingerprint != self.replay_fingerprint
            or right_token.dataset_fingerprint != self.replay_fingerprint
        ):
            raise ValueError("cross-dataset v0.23.3 transport is forbidden.")
        return left_token, right_token

    def snapshot_overlap(self, left, right):
        left_token, right_token = self._tokens(left, right)
        return self.dataset.overlaps[
            left_token.record_index, right_token.record_index
        ].copy()

    def snapshot_transport(self, left, right):
        left_token, right_token = self._tokens(left, right)
        return self.dataset.transport(
            left_token.record_index, right_token.record_index
        )

    def diagnostics_dict(self):
        overlap = self.dataset.overlap_diagnostics()
        return {
            "provider": "FileBackedMolecularSOCProviderV233",
            "dataset_fingerprint": self.replay_fingerprint,
            "legacy_dataset_fingerprint": self.dataset.legacy_dataset_fingerprint,
            "records": int(len(self.dataset.q)),
            "exact_record_only": True,
            "raw_overlap_contract": OVERLAP_CONTRACT_ID_V233,
            "transport_contract": TRANSPORT_CONTRACT_ID_V233,
            "overlap_quality_policy": self.dataset.overlap_policy.as_dict(),
            "overlap_diagnostics": overlap.as_dict(),
            "nac_convention_fingerprint": self.dataset.nac_convention.fingerprint(),
            "provider_numerical_identity_fingerprint": (
                self.numerical_identity_v233.fingerprint()
            ),
            "all_electronic_calculations_converged": bool(
                np.all(self.dataset.converged)
            ),
        }
