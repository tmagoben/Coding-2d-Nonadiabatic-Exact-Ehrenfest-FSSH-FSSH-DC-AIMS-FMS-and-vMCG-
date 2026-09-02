"""External molecular-SOC snapshot intake campaign for v0.24.0."""

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
from scipy.linalg import expm

from .external_soc_admission_v240 import (
    ExternalSOCAdmissionPolicyV240,
    audit_external_soc_snapshot_v240,
    require_external_soc_snapshot_v240,
)
from .external_soc_dynamics_v240 import (
    FrozenSnapshotCheckpointV240,
    preview_frozen_snapshot_dynamics_v240,
    run_admitted_external_soc_dynamics_v240,
)
from .external_soc_validation_v240 import (
    EXTERNAL_VALIDATION_SCHEMA_V240,
    audit_external_soc_validation_v240,
)
from .openmolcas_rassi_protocol_v240 import (
    OPENMOLCAS_ADAPTER_NAME_V240,
    OPENMOLCAS_ADAPTER_VERSION_V240,
    OPENMOLCAS_EXPORT_SCHEMA_V240,
    OPENMOLCAS_MANIFEST_SCHEMA_V240,
    openmolcas_protocol_from_dict_v240,
    water_rassi_so_protocol_v240,
)
from .openmolcas_rassi_snapshot_v240 import (
    OPENMOLCAS_EXPORT_NAME_V240,
    OPENMOLCAS_HDF5_NAME_V240,
    OPENMOLCAS_INPUT_NAME_V240,
    OPENMOLCAS_MANIFEST_NAME_V240,
    OPENMOLCAS_OUTPUT_NAME_V240,
    OPENMOLCAS_VALIDATION_NAME_V240,
    OPENMOLCAS_VALIDATION_ARTIFACT_DIRECTORY_V240,
    PROTOCOL_FIXTURE_MARKER_V240,
    OpenMolcasRASSISnapshotParserV240,
    sha256_file_v240,
)
from .soc_derivative_evidence_v240 import audit_external_soc_derivatives_v240
from .v233_benchmark import run_v0233_release_benchmark


@dataclass(frozen=True)
class V240AcceptanceThresholds:
    expected_inherited_gates: int = 208
    expected_new_gates: int = 48
    expected_total_gates: int = 256


def _json_bytes_v240(payload):
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _write_json_v240(path, payload):
    Path(path).write_bytes(_json_bytes_v240(payload))


def _digest_text_v240(value):
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()


def _raises_v240(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _fixture_components_v240(geometry, reference):
    displacement = np.asarray(geometry, dtype=float).reshape(-1) - np.asarray(
        reference, dtype=float
    ).reshape(-1)
    coordinate_weights = np.linspace(-0.7, 0.9, len(displacement))
    singlet_gradient = 2.0e-3 * coordinate_weights
    triplet_gradient = -1.5e-3 * coordinate_weights[::-1]
    singlet = -75.0 + np.dot(singlet_gradient, displacement) + 0.001 * np.dot(
        displacement, displacement
    )
    triplet = -74.8 + np.dot(triplet_gradient, displacement) + 0.0008 * np.dot(
        displacement, displacement
    )
    H0 = np.diag([singlet, triplet, triplet, triplet]).astype(complex)
    scalar = float(np.dot(np.linspace(0.02, -0.01, len(displacement)), displacement))
    lam = 0.0018 * ((0.73 + scalar) - 0.24j)
    mu = 0.0018 * (0.41 - 0.5 * scalar)
    Hso = np.zeros((4, 4), dtype=complex)
    Hso[0, 1:] = np.asarray([lam, 1j * mu, np.conj(lam)])
    Hso[1:, 0] = np.conj(Hso[0, 1:])
    return H0, Hso


def _record_geometry_v240(protocol, record_id):
    geometry = np.asarray(protocol.reference_geometry_bohr, dtype=float).copy()
    if record_id == "reference":
        return geometry
    pieces = record_id.split("_")
    coordinate = int(pieces[0][1:])
    step_index = int(pieces[1][1:])
    sign = -1.0 if pieces[2] == "minus" else 1.0
    geometry.reshape(-1)[coordinate] += sign * protocol.displacement_steps_bohr[step_index]
    return geometry


def build_v240_protocol_fixture(directory):
    """Create a deterministic, conspicuously non-ab-initio parser fixture."""

    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    protocol = water_rassi_so_protocol_v240()
    records = []
    output_text = (
        PROTOCOL_FIXTURE_MARKER_V240
        + "\nSynthetic parser exercise only; RASSI and SPIN-ORBIT words are not output.\n"
    )
    for record_id in protocol.expected_record_ids():
        record_directory = root / "records" / record_id
        record_directory.mkdir(parents=True)
        input_path = record_directory / OPENMOLCAS_INPUT_NAME_V240
        output_path = record_directory / OPENMOLCAS_OUTPUT_NAME_V240
        hdf5_path = record_directory / OPENMOLCAS_HDF5_NAME_V240
        export_path = record_directory / OPENMOLCAS_EXPORT_NAME_V240
        geometry = _record_geometry_v240(protocol, record_id)
        geometry_lines = "\n".join(
            f"{symbol} {row[0]:.15f} {row[1]:.15f} {row[2]:.15f}"
            for symbol, row in zip(protocol.atom_symbols, geometry)
        )
        input_text = f"""&GATEWAY
Coord
3
water displaced record
* GND_GEOMETRY_BOHR_BEGIN
{geometry_lines}
* GND_GEOMETRY_BOHR_END
Basis=ANO-RCC-VDZP
Unit=Bohr
Group=Nosym
Rela=R02O
&SEWARD
AMFI
&SCF
&RASSCF
FileOrb=$Project.ScfOrb
NACTEL=8 0 0
RAS2=6
Inactive=1
Spin=1
CIRoot=1 1 1
&CASPT2
>> COPY $Project.JobMix S.JobMix
&RASSCF
FileOrb=$Project.ScfOrb
NACTEL=8 0 0
RAS2=6
Inactive=1
Spin=3
CIRoot=1 1 1
&CASPT2
>> COPY $Project.JobMix T.JobMix
>> COPY S.JobMix JOB001
>> COPY T.JobMix JOB002
&RASSI
NROFJOBIPH=2 1 1;1;1
SPINORBIT
EJOB
"""
        input_path.write_text(input_text, encoding="utf-8")
        output_path.write_text(output_text, encoding="utf-8")
        hdf5_path.write_bytes(b"GND protocol fixture, deliberately not HDF5\n")
        H0, Hso = _fixture_components_v240(
            geometry, protocol.reference_geometry_bohr
        )
        export = {
            "schema": OPENMOLCAS_EXPORT_SCHEMA_V240,
            "record_id": record_id,
            "protocol_fingerprint": protocol.fingerprint(),
            "input_sha256": sha256_file_v240(input_path),
            "output_sha256": sha256_file_v240(output_path),
            "rassi_h5_sha256": sha256_file_v240(hdf5_path),
            "geometry_bohr": geometry.tolist(),
            "state_labels": list(protocol.state_order),
            "H_spin_free_real": H0.real.tolist(),
            "H_spin_free_imag": H0.imag.tolist(),
            "H_soc_real": Hso.real.tolist(),
            "H_soc_imag": Hso.imag.tolist(),
            "reference_overlap_real": np.eye(4).tolist(),
            "reference_overlap_imag": np.zeros((4, 4)).tolist(),
            "convergence": {
                "gateway": True,
                "seward": True,
                "scf": True,
                "rasscf_singlet": True,
                "rasscf_triplet": True,
                "caspt2_singlet": True,
                "caspt2_triplet": True,
                "rassi_so": True,
            },
        }
        _write_json_v240(export_path, export)
        records.append(
            {
                "record_id": record_id,
                "relative_directory": f"records/{record_id}",
                "input_sha256": sha256_file_v240(input_path),
                "output_sha256": sha256_file_v240(output_path),
                "rassi_h5_sha256": sha256_file_v240(hdf5_path),
                "export_sha256": sha256_file_v240(export_path),
            }
        )
    reference_Hso = _fixture_components_v240(
        protocol.reference_geometry_bohr, protocol.reference_geometry_bohr
    )[1]
    basis_values = np.asarray(
        [0.98 * reference_Hso, 0.99999 * reference_Hso, reference_Hso]
    )
    method_values = np.asarray(
        [0.97 * reference_Hso, 0.99998 * reference_Hso, reference_Hso]
    )
    eigenvalues = np.linalg.eigvalsh(reference_Hso)
    displaced_count = len(records) - 1
    validation = {
        "schema": EXTERNAL_VALIDATION_SCHEMA_V240,
        "source_kind": "protocol_fixture",
        "protocol_fingerprint": protocol.fingerprint(),
        "reference_export_sha256": records[0]["export_sha256"],
        "independent_backend_name": "independent analytic protocol fixture",
        "independent_backend_version": "1",
        "independent_artifact_sha256": _digest_text_v240("fixture-reference"),
        "state_order": list(protocol.state_order),
        "reference_soc_real": reference_Hso.real.tolist(),
        "reference_soc_imag": reference_Hso.imag.tolist(),
        "independent_soc_real": (reference_Hso * (1.0 + 1.0e-7)).real.tolist(),
        "independent_soc_imag": (reference_Hso * (1.0 + 1.0e-7)).imag.tolist(),
        "reference_tolerance_hartree": 1.0e-8,
        "basis_labels": ["ANO-RCC-MB", "ANO-RCC-VDZP", "ANO-RCC-VTZP"],
        "basis_soc_real": basis_values.real.tolist(),
        "basis_soc_imag": basis_values.imag.tolist(),
        "basis_artifact_sha256": [
            _digest_text_v240(f"basis-{index}") for index in range(3)
        ],
        "basis_tolerance_hartree": 1.0e-7,
        "method_labels": ["CASSCF", "MS-CASPT2", "SS-CASPT2-EJOB"],
        "method_soc_real": method_values.real.tolist(),
        "method_soc_imag": method_values.imag.tolist(),
        "method_artifact_sha256": [
            _digest_text_v240(f"method-{index}") for index in range(3)
        ],
        "method_tolerance_hartree": 1.0e-7,
        "frame_base_eigenvalues_hartree": eigenvalues.tolist(),
        "frame_translated_eigenvalues_hartree": eigenvalues.tolist(),
        "frame_rotated_eigenvalues_hartree": eigenvalues.tolist(),
        "frame_artifact_sha256": [
            _digest_text_v240(f"frame-{index}") for index in range(3)
        ],
        "frame_tolerance_hartree": 1.0e-9,
        "tracking_minimum_singular_values": [1.0] * displaced_count,
        "tracking_maximum_competing_leakage": [0.0] * displaced_count,
        "tracking_assignment_margins": [1.0] * displaced_count,
        "tracking_artifact_sha256": [
            _digest_text_v240(f"tracking-{index}")
            for index in range(displaced_count)
        ],
    }
    validation_path = root / OPENMOLCAS_VALIDATION_NAME_V240
    _write_json_v240(validation_path, validation)
    validation_artifact_directory = (
        root / OPENMOLCAS_VALIDATION_ARTIFACT_DIRECTORY_V240
    )
    validation_artifact_directory.mkdir()
    validation_artifact_payloads = [
        "fixture-reference",
        *[f"basis-{index}" for index in range(3)],
        *[f"method-{index}" for index in range(3)],
        *[f"frame-{index}" for index in range(3)],
        *[f"tracking-{index}" for index in range(displaced_count)],
    ]
    for payload in validation_artifact_payloads:
        digest = _digest_text_v240(payload)
        (validation_artifact_directory / f"{digest}.artifact").write_bytes(
            payload.encode("utf-8")
        )
    environment_sha256 = _digest_text_v240(
        "protocol fixture: no OpenMolcas runtime environment"
    )
    manifest = {
        "schema": OPENMOLCAS_MANIFEST_SCHEMA_V240,
        "source_kind": "protocol_fixture",
        "protocol": protocol.as_dict(),
        "adapter_name": OPENMOLCAS_ADAPTER_NAME_V240,
        "adapter_version": OPENMOLCAS_ADAPTER_VERSION_V240,
        "exporter_name": "gnd-rassi-hdf5-exporter",
        "exporter_version": "0.24.0",
        "environment_sha256": environment_sha256,
        "validation_artifact": OPENMOLCAS_VALIDATION_NAME_V240,
        "validation_sha256": sha256_file_v240(validation_path),
        "records": records,
    }
    manifest_path = root / OPENMOLCAS_MANIFEST_NAME_V240
    _write_json_v240(manifest_path, manifest)
    return root


def _fixture_policy_v240(directory, protocol=None):
    directory = Path(directory)
    protocol = water_rassi_so_protocol_v240() if protocol is None else protocol
    manifest = json.loads(
        (directory / OPENMOLCAS_MANIFEST_NAME_V240).read_text(encoding="utf-8")
    )
    return ExternalSOCAdmissionPolicyV240(
        expected_protocol=protocol,
        expected_soc_convention=protocol.soc_convention(),
        trusted_parser_type=OpenMolcasRASSISnapshotParserV240,
        expected_manifest_sha256=sha256_file_v240(
            directory / OPENMOLCAS_MANIFEST_NAME_V240
        ),
        expected_environment_sha256=manifest["environment_sha256"],
        expected_exporter_name=manifest["exporter_name"],
        expected_exporter_version=manifest["exporter_version"],
    ).validate()


def _protocol_controls_v240():
    protocol = water_rassi_so_protocol_v240()
    restored = openmolcas_protocol_from_dict_v240(protocol.as_dict())
    convention = protocol.soc_convention()
    model_space = protocol.model_space()
    return {
        "openmolcas_version_pinned": protocol.backend_version == "26.06",
        "water_nuclear_identity_pinned": bool(
            protocol.atom_symbols == ("O", "H", "H")
            and len(protocol.isotope_masses_amu) == 3
        ),
        "neutral_ten_electron_sector": bool(
            protocol.charge == 0 and protocol.electron_count == 10
        ),
        "complete_singlet_triplet_space": bool(
            model_space.complete_multiplets and model_space.nstate == 4
        ),
        "exact_spin_component_order": tuple(convention.state_order)
        == protocol.state_order,
        "nine_cartesian_coordinates": protocol.coordinate_dimension == 9,
        "three_decreasing_displacements": protocol.displacement_steps_bohr
        == (0.004, 0.002, 0.001),
        "amfi_soc_operator_frozen": "AMFI" in protocol.soc_operator,
        "dkh2_scalar_method_frozen": "DKH2" in protocol.scalar_relativistic_method,
        "caspt2_ejob_energy_source_frozen": bool(
            "CASPT2" in protocol.dynamic_correlation_method
            and protocol.rassi_energy_source == "EJOB"
        ),
        "cross_geometry_overlap_not_mislabeled_rassi": "not inferred" in protocol.cross_geometry_overlap_method,
        "protocol_roundtrip_and_fingerprint": bool(
            restored == protocol and restored.fingerprint() == protocol.fingerprint()
        ),
        "protocol_fingerprint": protocol.fingerprint(),
    }


def _parser_derivative_controls_v240(directory):
    parser = OpenMolcasRASSISnapshotParserV240()
    bundle = parser.parse_bundle_v240(directory)
    derivative = audit_external_soc_derivatives_v240(bundle)
    validation = audit_external_soc_validation_v240(bundle)
    return {
        "strict_manifest_parsed": bundle.manifest.schema == OPENMOLCAS_MANIFEST_SCHEMA_V240,
        "complete_55_record_inventory": len(bundle.records) == 55,
        "exact_artifact_inventory": bundle.exact_artifact_inventory,
        "typed_parser_executed": bundle.parser_executed,
        "fixture_source_preserved": bundle.source_kind == "protocol_fixture",
        "fixture_not_native_execution": not bundle.native_openmolcas_execution,
        "derivative_evidence_passes": derivative.passed,
        "complete_centered_geometry_pairs": derivative.checks["complete_centered_geometry_pairs"],
        "all_component_calculations_converged": derivative.checks["all_calculations_converged"],
        "component_and_derivative_hermiticity": bool(
            derivative.checks["component_hermiticity"]
            and derivative.checks["derivative_hermiticity"]
        ),
        "finite_manifold_polar_transport": derivative.checks["finite_manifold_transport"],
        "complete_manifold_tracking": derivative.checks["complete_manifold_tracking"],
        "independent_reference_agreement": validation.checks["independent_reference_agreement"],
        "basis_convergence": validation.checks["basis_convergence"],
        "method_convergence": validation.checks["method_convergence"],
        "frame_and_tracking_validation": bool(
            validation.checks["translation_invariance"]
            and validation.checks["rotation_invariance"]
            and validation.checks["tracking_retention"]
            and validation.checks["tracking_leakage"]
            and validation.checks["tracking_margin"]
        ),
        "bundle": bundle,
        "derivative_report": derivative.as_dict(),
        "validation_report": validation.as_dict(),
    }


def _admission_controls_v240(directory):
    directory = Path(directory)
    parser = OpenMolcasRASSISnapshotParserV240()
    policy = _fixture_policy_v240(directory)
    audit = audit_external_soc_snapshot_v240(directory, policy=policy, parser=parser)

    relabel = directory.parent / "relabel"
    shutil.copytree(directory, relabel)
    relabel_manifest_path = relabel / OPENMOLCAS_MANIFEST_NAME_V240
    relabel_manifest = json.loads(relabel_manifest_path.read_text(encoding="utf-8"))
    relabel_manifest["source_kind"] = "external_ab_initio_snapshot"
    _write_json_v240(relabel_manifest_path, relabel_manifest)

    corrupted_output = directory.parent / "corrupted-output"
    shutil.copytree(directory, corrupted_output)
    with (corrupted_output / "records" / "reference" / OPENMOLCAS_OUTPUT_NAME_V240).open(
        "a", encoding="utf-8"
    ) as handle:
        handle.write("corruption\n")

    geometry_mismatch = directory.parent / "geometry-mismatch"
    shutil.copytree(directory, geometry_mismatch)
    mismatch_input = geometry_mismatch / "records" / "reference" / OPENMOLCAS_INPUT_NAME_V240
    mismatch_text = mismatch_input.read_text(encoding="utf-8").replace(
        "O 0.000000000000000 0.000000000000000 0.000000000000000",
        "O 0.100000000000000 0.000000000000000 0.000000000000000",
        1,
    )
    mismatch_input.write_text(mismatch_text, encoding="utf-8")
    mismatch_export_path = (
        geometry_mismatch / "records" / "reference" / OPENMOLCAS_EXPORT_NAME_V240
    )
    mismatch_export = json.loads(mismatch_export_path.read_text(encoding="utf-8"))
    mismatch_export["input_sha256"] = sha256_file_v240(mismatch_input)
    _write_json_v240(mismatch_export_path, mismatch_export)
    mismatch_manifest_path = geometry_mismatch / OPENMOLCAS_MANIFEST_NAME_V240
    mismatch_manifest = json.loads(mismatch_manifest_path.read_text(encoding="utf-8"))
    mismatch_manifest["records"][0]["input_sha256"] = sha256_file_v240(mismatch_input)
    mismatch_manifest["records"][0]["export_sha256"] = sha256_file_v240(
        mismatch_export_path
    )
    _write_json_v240(mismatch_manifest_path, mismatch_manifest)

    unknown = directory.parent / "unknown-artifact"
    shutil.copytree(directory, unknown)
    (unknown / "unexpected.bin").write_bytes(b"unexpected")

    bad_validation = directory.parent / "bad-validation"
    shutil.copytree(directory, bad_validation)
    with (bad_validation / OPENMOLCAS_VALIDATION_NAME_V240).open("a", encoding="utf-8") as handle:
        handle.write(" ")

    bad_validation_blob = directory.parent / "bad-validation-blob"
    shutil.copytree(directory, bad_validation_blob)
    first_blob = next(
        iter(
            sorted(
                (
                    bad_validation_blob
                    / OPENMOLCAS_VALIDATION_ARTIFACT_DIRECTORY_V240
                ).glob("*.artifact")
            )
        )
    )
    with first_blob.open("ab") as handle:
        handle.write(b"corruption")

    wrong_environment_policy = replace(policy, expected_environment_sha256="0" * 64)
    wrong_environment = audit_external_soc_snapshot_v240(
        directory, policy=wrong_environment_policy, parser=parser
    )
    wrong_convention = replace(
        policy.expected_soc_convention,
        prefactor_convention="untrusted adapter-side factor of two",
    )
    wrong_convention_policy = replace(policy, expected_soc_convention=wrong_convention)

    class ParserSubclass(OpenMolcasRASSISnapshotParserV240):
        pass

    subclass_audit = audit_external_soc_snapshot_v240(
        directory, policy=policy, parser=ParserSubclass()
    )
    return {
        "protocol_fixture_passes_protocol_audit": audit.protocol_passed,
        "protocol_fixture_not_external_admitted": not audit.external_snapshot_admitted,
        "live_backend_remains_closed": not audit.live_backend_admitted,
        "production_requirement_rejects_fixture": _raises_v240(
            lambda: require_external_soc_snapshot_v240(
                directory, policy=policy, parser=parser
            ),
            (ValueError,),
            "not admitted",
        ),
        "synthetic_relabel_rejected": _raises_v240(
            lambda: parser.parse_bundle_v240(relabel),
            (ValueError,),
            "cannot be relabeled",
        ),
        "native_output_and_input_geometry_corruption_rejected": bool(
            _raises_v240(
                lambda: parser.parse_bundle_v240(corrupted_output),
                (ValueError,),
                "digest mismatch",
            )
            and _raises_v240(
                lambda: parser.parse_bundle_v240(geometry_mismatch),
                (ValueError,),
                "not bound to the native input",
            )
        ),
        "unknown_artifact_rejected": _raises_v240(
            lambda: parser.parse_bundle_v240(unknown),
            (ValueError,),
            "unknown or missing",
        ),
        "environment_trust_anchor_enforced": not wrong_environment.protocol_passed,
        "soc_convention_trust_anchor_enforced": _raises_v240(
            wrong_convention_policy.validate,
            (ValueError,),
            "disagree",
        ),
        "parser_subclass_rejected": bool(
            not subclass_audit.checks["trusted_parser_exact_type"]
            and not subclass_audit.protocol_passed
        ),
        "validation_json_and_raw_blob_corruption_rejected": bool(
            _raises_v240(
                lambda: parser.parse_bundle_v240(bad_validation),
                (ValueError,),
                "artifact digest mismatch",
            )
            and _raises_v240(
                lambda: parser.parse_bundle_v240(bad_validation_blob),
                (ValueError,),
                "raw artifact digest mismatch",
            )
        ),
        "external_claim_and_native_crosscheck_are_explicitly_false": bool(
            audit.source_kind == "protocol_fixture"
            and not audit.checks["native_numeric_crosscheck"]
            and not audit.external_snapshot_admitted
        ),
        "audit": audit.as_dict(),
    }


def _dynamics_controls_v240(bundle, admission_audit):
    initial = np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex)
    dt = 0.25
    full = preview_frozen_snapshot_dynamics_v240(
        bundle, initial, time_step_au=dt, steps=12
    )
    spin_free = preview_frozen_snapshot_dynamics_v240(
        bundle, initial, time_step_au=dt, steps=12, soc_enabled=False
    )
    reference = bundle.record_map["reference"]
    expected = expm(-1j * 12 * dt * reference.H_total) @ initial
    expected_spin_free = expm(-1j * 12 * dt * reference.H_spin_free) @ initial
    first = preview_frozen_snapshot_dynamics_v240(
        bundle, initial, time_step_au=dt, steps=5
    )
    second = preview_frozen_snapshot_dynamics_v240(
        bundle,
        initial,
        time_step_au=dt,
        steps=7,
        checkpoint=first.checkpoint,
    )
    mismatched_checkpoint = FrozenSnapshotCheckpointV240(
        bundle_fingerprint="0" * 64,
        step=first.checkpoint.step,
        time_au=first.checkpoint.time_au,
        coefficients=first.checkpoint.coefficients,
    )
    return {
        "preview_preserves_fixture_label": full.evidence_class == "protocol_fixture",
        "unitary_norm_conservation": float(np.max(np.abs(full.norms - 1.0))) < 1.0e-12,
        "static_matrix_exponential_reference": np.allclose(full.coefficients[-1], expected, atol=1.0e-12),
        "zero_soc_mode_is_distinct": not np.allclose(full.coefficients[-1], spin_free.coefficients[-1]),
        "zero_soc_matrix_exponential_reference": np.allclose(
            spin_free.coefficients[-1], expected_spin_free, atol=1.0e-12
        ),
        "deterministic_checkpoint_restart": np.array_equal(
            second.coefficients[-1], full.coefficients[-1]
        ),
        "production_dynamics_rejects_fixture": _raises_v240(
            lambda: run_admitted_external_soc_dynamics_v240(
                bundle,
                admission_audit,
                initial,
                time_step_au=dt,
                steps=1,
            ),
            (ValueError,),
            "admitted external snapshot",
        ),
        "checkpoint_bundle_identity_enforced": _raises_v240(
            lambda: preview_frozen_snapshot_dynamics_v240(
                bundle,
                initial,
                time_step_au=dt,
                steps=1,
                checkpoint=mismatched_checkpoint,
            ),
            (ValueError,),
            "identities differ",
        ),
        "norm_error": float(np.max(np.abs(full.norms - 1.0))),
        "restart_error": float(np.linalg.norm(second.coefficients[-1] - full.coefficients[-1])),
    }


def run_v0240_release_benchmark(thresholds=V240AcceptanceThresholds()):
    inherited = run_v0233_release_benchmark()
    protocol = _protocol_controls_v240()
    with tempfile.TemporaryDirectory(prefix="gnd-v240-release-") as temporary:
        temporary = Path(temporary)
        fixture = build_v240_protocol_fixture(temporary / "fixture")
        parser_derivative = _parser_derivative_controls_v240(fixture)
        admission = _admission_controls_v240(fixture)
        dynamics = _dynamics_controls_v240(
            parser_derivative["bundle"],
            audit_external_soc_snapshot_v240(
                fixture,
                policy=_fixture_policy_v240(fixture),
                parser=OpenMolcasRASSISnapshotParserV240(),
            ),
        )
    new_checks = {
        **{
            f"protocol::{name}": bool(value)
            for name, value in protocol.items()
            if name != "protocol_fingerprint"
        },
        **{
            f"parser_derivative::{name}": bool(value)
            for name, value in parser_derivative.items()
            if name not in {"bundle", "derivative_report", "validation_report"}
        },
        **{
            f"admission::{name}": bool(value)
            for name, value in admission.items()
            if name != "audit"
        },
        **{
            f"dynamics::{name}": bool(value)
            for name, value in dynamics.items()
            if name not in {"norm_error", "restart_error"}
        },
    }
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError(
            f"v0.24.0 defines {len(new_checks)} rather than exactly "
            f"{thresholds.expected_new_gates} new gates."
        )
    inherited_checks = {
        f"inherited_v0233::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.24.0 must inherit exactly 208 v0.23.3 gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.24.0 must define exactly 256 cumulative gates.")
    return {
        "release": "v0.24.0",
        "theme": "fail-closed OpenMolcas RASSI-SO external snapshot intake",
        "protocol_controls": protocol,
        "parser_and_derivative_controls": {
            name: value
            for name, value in parser_derivative.items()
            if name != "bundle"
        },
        "admission_controls": admission,
        "dynamics_controls": dynamics,
        "claims": {
            "openmolcas_rassi_so_protocol_frozen": True,
            "strict_bundle_artifact_parser_validated": True,
            "transported_cartesian_soc_derivative_protocol_validated": True,
            "independent_accuracy_evidence_schema_validated": True,
            "admission_bound_frozen_snapshot_dynamics_validated": True,
            "protocol_fixture_validated": True,
            "external_molecular_SOC_snapshot_admitted": False,
            "live_molecular_SOC_backend_admitted": False,
            "ab_initio_SOC_validated": False,
            "openmolcas_runtime_executed": False,
            "native_openmolcas_numeric_crosscheck_implemented": False,
        },
        "inherited_v0233": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
