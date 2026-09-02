"""Raw-evidence and runtime-attestation admission campaign for v0.23.1."""

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import shutil
import tempfile
import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
)
from .molecular_soc_admission_v230 import audit_molecular_soc_provider_v230
from .molecular_soc_admission_v231 import audit_molecular_soc_provider_v231
from .molecular_soc_contract_v230 import (
    MolecularSOCAdmissionContractV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCCapabilitiesV230,
    MolecularSOCValidationEvidenceV230,
)
from .molecular_soc_dossier_v231 import (
    DOSSIER_NAME_V231,
    CalculationReceiptV231,
    MolecularSOCAdmissionDossierV231,
    RawArtifactRecordV231,
    load_molecular_soc_dossier_v231,
    write_molecular_soc_dossier_v231,
    write_raw_json_artifact_v231,
)
from .molecular_soc_evidence_v231 import (
    ConvergenceLadderObservationV231,
    DerivedEvidenceBundleV231,
    FrameInvarianceObservationV231,
    IndependentReferenceObservationV231,
    TrackingSpecificationV231,
)
from .molecular_soc_replay_v230 import (
    FileBackedMolecularSOCProviderV230,
    capture_molecular_soc_replay_v230,
)
from .pyscf_soc_adapter_v231 import (
    probe_pyscf_soc_adapter_v231,
    require_pyscf_soc_adapter_v231,
)
from .v230_benchmark import run_v0230_release_benchmark, v230_reference_coordinates


@dataclass(frozen=True)
class V231AcceptanceThresholds:
    expected_inherited_gates: int = 93
    expected_new_gates: int = 30
    expected_total_gates: int = 123


def _trajectory_capabilities_v231():
    return MolecularSOCCapabilitiesV230(
        static_soc=True,
        spin_free_derivatives=True,
        soc_derivatives=True,
        derivative_connections=True,
        cross_geometry_overlaps=True,
        deterministic_replay=True,
        analytic_soc_derivatives=True,
    ).validate()


def _scalar_soc_observable_v231(source, q):
    return float(np.linalg.norm(source.components(np.asarray(q, dtype=float)).H_soc))


def _receipt_artifacts_v231(
    directory,
    *,
    identity,
    record_id,
    role,
    q,
    output_payload,
    basis=None,
    method=None,
    converged=True,
    overwrite=False,
):
    input_name = f"{record_id}.input"
    output_name = f"{record_id}.output"
    input_record = write_raw_json_artifact_v231(
        directory,
        name=input_name,
        relative_path=f"raw/receipts/{record_id}.input.json",
        role="calculation_input",
        payload={
            "record_id": record_id,
            "role": role,
            "q_bohr": list(np.asarray(q, dtype=float)),
            "electronic_method": identity.electronic_method if method is None else method,
            "basis": identity.basis if basis is None else basis,
            "soc_operator": identity.soc_operator,
            "derivative_method": identity.derivative_method,
        },
        overwrite=overwrite,
    )
    output_record = write_raw_json_artifact_v231(
        directory,
        name=output_name,
        relative_path=f"raw/receipts/{record_id}.output.json",
        role="calculation_output",
        payload={
            "record_id": record_id,
            "role": role,
            "q_bohr": list(np.asarray(q, dtype=float)),
            "result": output_payload,
            "convergence": {
                "scf": bool(converged),
                "correlated": bool(converged),
                "soc": bool(converged),
                "derivatives": bool(converged),
                "overlaps": bool(converged),
            },
        },
        overwrite=overwrite,
    )
    receipt = CalculationReceiptV231(
        record_id=record_id,
        role=role,
        q_bohr=tuple(float(value) for value in np.asarray(q, dtype=float)),
        backend_name=identity.backend_name,
        backend_version=identity.backend_version,
        source_kind=identity.source_kind,
        electronic_method=identity.electronic_method if method is None else method,
        basis=identity.basis if basis is None else basis,
        soc_operator=identity.soc_operator,
        derivative_method=identity.derivative_method,
        input_artifact=input_name,
        output_artifact=output_name,
        scf_converged=bool(converged),
        correlated_converged=bool(converged),
        soc_converged=bool(converged),
        derivatives_converged=bool(converged),
        overlaps_converged=bool(converged),
    ).validate()
    return (input_record, output_record), receipt


def _tracking_specification_v231(odd, nrecord):
    if odd:
        labels = ("doublet_1", "doublet_2")
        groups = ((0, 1), (2, 3))
    else:
        labels = ("singlet", "triplet")
        groups = ((0,), (1, 2, 3))
    return TrackingSpecificationV231(
        manifold_labels=labels,
        manifold_state_indices=groups,
        record_edges=tuple((index, index + 1) for index in range(nrecord - 1)),
        overlap_threshold=0.80,
        margin_threshold=0.10,
    ).validate(nrecord=nrecord, nstate=4)


def build_v231_admission_bundle(
    directory,
    *,
    odd=False,
    source_kind="validation_fixture",
    overwrite=False,
):
    """Build a molecular-format fixture bundle without granting it a real claim."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = AnalyticDoubletSOCProviderV220() if odd else AnalyticSingletTripletSOCProviderV220()
    center = -0.11 if odd else 0.17
    coordinates = v230_reference_coordinates(center)
    sector = "Kramers-doublet" if odd else "singlet-triplet"
    backend_name = f"deterministic analytic {sector} raw-evidence fixture"
    backend_version = "1"
    electronic_method = f"analytic {sector} reference"
    basis_name = "closed-form model basis"
    soc_operator = "analytic time-reversal-even SOC"
    derivative_method = "analytic component derivatives"

    template = write_raw_json_artifact_v231(
        directory,
        name="calculation_template",
        relative_path="raw/calculation_template.json",
        role="calculation_template",
        payload={
            "backend": backend_name,
            "backend_version": backend_version,
            "source_kind": source_kind,
            "sector": sector,
            "coordinate_unit": "bohr",
            "energy_unit": "hartree",
        },
        overwrite=overwrite,
    )
    environment = write_raw_json_artifact_v231(
        directory,
        name="environment_lock",
        relative_path="raw/environment_lock.json",
        role="environment_lock",
        payload={
            "environment": "deterministic analytic validation fixture",
            "numpy_required": True,
            "real_electronic_structure_runtime": False,
        },
        overwrite=overwrite,
    )
    external_reference = write_raw_json_artifact_v231(
        directory,
        name="independent_reference",
        relative_path="raw/independent_reference.json",
        role="independent_reference",
        payload={
            "reference_id": "independent analytic negative-control table",
            "physical_molecular_reference": False,
            "sector": sector,
        },
        overwrite=overwrite,
    )
    artifacts = [template, environment, external_reference]

    real_source = source_kind in {"external_ab_initio_snapshot", "live_ab_initio"}
    identity = MolecularSOCBackendIdentityV230(
        backend_name=backend_name,
        backend_version=backend_version,
        source_kind=source_kind,
        electronic_method=electronic_method,
        basis=basis_name,
        charge=0,
        electron_count=3 if odd else 2,
        soc_operator=soc_operator,
        scalar_relativistic_method="none",
        derivative_method=derivative_method,
        active_space="complete four-state validation space",
        molecule_name=(
            "synthetically relabelled negative-control pseudo-molecule"
            if real_source
            else "analytic validation fixture"
        ),
        atom_symbols=(("H", "H") if real_source else ()),
        isotope_masses_amu=((1.007825, 1.007825) if real_source else ()),
        reference_geometry_bohr=(
            ((0.0, 0.0, -0.7), (0.0, 0.0, 0.7)) if real_source else ()
        ),
        calculation_input_sha256=template.sha256,
        environment_sha256=environment.sha256,
        extra={
            "molecular_accuracy_claim": False,
            "synthetic_relabel_negative_control": bool(real_source),
        },
    ).validate()

    receipts = []
    trajectory_ids = []
    for index, q in enumerate(coordinates):
        record_id = f"trajectory_{index:03d}"
        components = source.components(q).validate()
        records, receipt = _receipt_artifacts_v231(
            directory,
            identity=identity,
            record_id=record_id,
            role="trajectory",
            q=q,
            output_payload={
                "H_spin_free": components.H_spin_free,
                "K_spin_free": components.K_spin_free,
                "H_soc": components.H_soc,
                "K_soc": components.K_soc,
            },
            overwrite=overwrite,
        )
        artifacts.extend(records)
        receipts.append(receipt)
        trajectory_ids.append(record_id)

    scalar = _scalar_soc_observable_v231(source, [center])
    basis_labels = ("fixture-basis-small", "fixture-basis-medium", "fixture-basis-large")
    basis_values = ((scalar + 2.0e-4,), (scalar + 5.0e-6,), (scalar,))
    basis_artifacts = []
    for index, (label, value) in enumerate(zip(basis_labels, basis_values)):
        records, receipt = _receipt_artifacts_v231(
            directory,
            identity=identity,
            record_id=f"basis_{index:03d}",
            role="basis",
            q=[center],
            basis=label,
            output_payload={"soc_norm_hartree": value[0]},
            overwrite=overwrite,
        )
        artifacts.extend(records)
        receipts.append(receipt)
        basis_artifacts.append(receipt.output_artifact)

    method_labels = ("fixture-method-a", "fixture-method-b")
    method_values = ((scalar + 4.0e-6,), (scalar,))
    method_artifacts = []
    for index, (label, value) in enumerate(zip(method_labels, method_values)):
        records, receipt = _receipt_artifacts_v231(
            directory,
            identity=identity,
            record_id=f"method_{index:03d}",
            role="method",
            q=[center],
            method=label,
            output_payload={"soc_norm_hartree": value[0]},
            overwrite=overwrite,
        )
        artifacts.extend(records)
        receipts.append(receipt)
        method_artifacts.append(receipt.output_artifact)

    frame_values = ((scalar,), (scalar + 1.0e-9,), (scalar + 2.0e-9,))
    frame_roles = ("frame_base", "frame_translation", "frame_rotation")
    frame_artifacts = []
    for index, (role, value) in enumerate(zip(frame_roles, frame_values)):
        records, receipt = _receipt_artifacts_v231(
            directory,
            identity=identity,
            record_id=role,
            role=role,
            q=[center],
            output_payload={"soc_norm_hartree": value[0]},
            overwrite=overwrite,
        )
        artifacts.extend(records)
        receipts.append(receipt)
        frame_artifacts.append(receipt.output_artifact)

    reference = IndependentReferenceObservationV231(
        reference_id="independent analytic negative-control table",
        observable="SOC Frobenius norm",
        unit="hartree",
        value_shape=(1,),
        computed_values=(scalar,),
        reference_values=(scalar + 1.0e-6,),
        computed_artifact=frame_artifacts[0],
        reference_artifact=external_reference.name,
        tolerance=1.0e-5,
    ).validate()
    basis = ConvergenceLadderObservationV231(
        kind="basis",
        labels=basis_labels,
        observable="SOC Frobenius norm",
        unit="hartree",
        value_shape=(1,),
        values=basis_values,
        source_artifacts=tuple(basis_artifacts),
        tolerance=1.0e-5,
    ).validate()
    method = ConvergenceLadderObservationV231(
        kind="method",
        labels=method_labels,
        observable="SOC Frobenius norm",
        unit="hartree",
        value_shape=(1,),
        values=method_values,
        source_artifacts=tuple(method_artifacts),
        tolerance=1.0e-5,
    ).validate()
    angle = 0.37
    rotation = (
        (float(np.cos(angle)), float(-np.sin(angle)), 0.0),
        (float(np.sin(angle)), float(np.cos(angle)), 0.0),
        (0.0, 0.0, 1.0),
    )
    frame = FrameInvarianceObservationV231(
        observable="SOC Frobenius norm",
        unit="hartree",
        value_shape=(1,),
        base_values=frame_values[0],
        translated_values=frame_values[1],
        rotated_values=frame_values[2],
        expected_rotated_values=frame_values[0],
        translation_bohr=(0.3, -0.2, 0.1),
        rotation_matrix=rotation,
        source_artifacts=tuple(frame_artifacts),
        tolerance=1.0e-8,
    ).validate()
    tracking = _tracking_specification_v231(odd, len(coordinates))
    raw_evidence = DerivedEvidenceBundleV231(
        reference=reference,
        basis=basis,
        method=method,
        frame=frame,
        tracking=tracking,
    ).validate()
    predicted_tracking = {
        "minimum_overlap": 1.0,
        "minimum_margin": 1.0,
    }
    evidence_v230 = raw_evidence.reference
    summarized_evidence = MolecularSOCValidationEvidenceV230(
        independent_reference_id=evidence_v230.reference_id,
        independent_reference_error=evidence_v230.error,
        independent_reference_tolerance=evidence_v230.tolerance,
        basis_levels=basis.labels,
        basis_changes=basis.changes,
        basis_tolerance=basis.tolerance,
        method_levels=method.labels,
        method_changes=method.changes,
        method_tolerance=method.tolerance,
        translation_residual=frame.translation_residual,
        rotation_residual=frame.rotation_residual,
        frame_invariance_tolerance=frame.tolerance,
        tracking_minimum_overlap=predicted_tracking["minimum_overlap"],
        tracking_minimum_margin=predicted_tracking["minimum_margin"],
        tracking_overlap_threshold=tracking.overlap_threshold,
        tracking_margin_threshold=tracking.margin_threshold,
    ).validate()
    contract = MolecularSOCAdmissionContractV230(
        capabilities=_trajectory_capabilities_v231(),
        identity=identity,
        evidence=summarized_evidence,
        state_tracking_policy="connected physical-manifold overlap graph",
        coordinate_definition="one explicit generalized coordinate in bohr",
        all_electronic_calculations_converged=True,
    ).validate(source.soc_symmetry_contract)
    replay_directory = directory / "replay"
    convergence_flags = np.ones(len(coordinates), dtype=bool) if real_source else None
    dataset = capture_molecular_soc_replay_v230(
        replay_directory,
        source,
        coordinates,
        contract,
        convergence_flags=convergence_flags,
        overwrite=overwrite,
    )
    dossier = MolecularSOCAdmissionDossierV231(
        replay_dataset_fingerprint=dataset.dataset_fingerprint,
        calculation_template_artifact=template.name,
        environment_artifact=environment.name,
        artifacts=tuple(artifacts),
        receipts=tuple(receipts),
        trajectory_record_ids=tuple(trajectory_ids),
        evidence=raw_evidence,
        runtime_attestation=None,
    )
    derived = dossier.derived_v230_evidence(dataset)
    if derived.as_dict() != contract.evidence.as_dict():
        raise AssertionError("fixture raw evidence does not reproduce its v0.23.0 summary.")
    dossier = write_molecular_soc_dossier_v231(
        directory,
        dossier,
        dataset=dataset,
        identity=identity,
        overwrite=overwrite,
    )
    provider = FileBackedMolecularSOCProviderV230(replay_directory)
    return {
        "directory": directory,
        "dataset": dataset,
        "provider": provider,
        "dossier": dossier,
        "dossier_path": directory / DOSSIER_NAME_V231,
        "center": np.asarray([center]),
    }


def _expect_value_error_v231(callable_object, text=None):
    try:
        callable_object()
    except (ValueError, FileNotFoundError) as exc:
        return text is None or text in str(exc)
    return False


def _campaign_controls_v231(work_directory):
    work_directory = Path(work_directory)
    even_first = build_v231_admission_bundle(work_directory / "even_first")
    even_second = build_v231_admission_bundle(work_directory / "even_second")
    odd = build_v231_admission_bundle(work_directory / "odd", odd=True)
    even_report = audit_molecular_soc_provider_v231(
        even_first["provider"],
        even_first["center"],
        even_first["dossier_path"],
        requirement="protocol",
    )
    odd_report = audit_molecular_soc_provider_v231(
        odd["provider"], odd["center"], odd["dossier_path"], requirement="protocol"
    )
    even_real = audit_molecular_soc_provider_v231(
        even_first["provider"],
        even_first["center"],
        even_first["dossier_path"],
        requirement="real",
    )
    odd_real = audit_molecular_soc_provider_v231(
        odd["provider"], odd["center"], odd["dossier_path"], requirement="real"
    )
    first_inventory = {
        record.name: record.sha256 for record in even_first["dossier"].artifacts
    }
    second_inventory = {
        record.name: record.sha256 for record in even_second["dossier"].artifacts
    }

    corrupted_artifact = work_directory / "corrupted_artifact"
    shutil.copytree(even_first["directory"], corrupted_artifact)
    target = corrupted_artifact / even_first["dossier"].artifacts[-1].relative_path
    target.write_bytes(target.read_bytes() + b"corruption")
    artifact_corruption_rejected = _expect_value_error_v231(
        lambda: load_molecular_soc_dossier_v231(
            corrupted_artifact,
            dataset=even_first["dataset"],
            identity=even_first["provider"].molecular_soc_contract.identity,
        )
    )

    corrupted_dossier = work_directory / "corrupted_dossier"
    shutil.copytree(even_first["directory"], corrupted_dossier)
    dossier_path = corrupted_dossier / DOSSIER_NAME_V231
    payload = json.loads(dossier_path.read_text(encoding="utf-8"))
    payload["trajectory_record_ids"][0] = "changed_without_fingerprint"
    dossier_path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    dossier_corruption_rejected = _expect_value_error_v231(
        lambda: load_molecular_soc_dossier_v231(dossier_path)
    )

    base_dossier = even_first["dossier"]
    missing_receipt = replace(base_dossier, receipts=base_dossier.receipts[1:])
    missing_receipt_rejected = _expect_value_error_v231(
        lambda: missing_receipt.validate(
            bundle_directory=even_first["directory"],
            dataset=even_first["dataset"],
            identity=even_first["provider"].molecular_soc_contract.identity,
        ),
        "trajectory_record_ids",
    )

    unconverged_directory = work_directory / "unconverged"
    shutil.copytree(even_first["directory"], unconverged_directory)
    unconverged_receipts = list(base_dossier.receipts)
    unconverged_receipts[0] = replace(
        unconverged_receipts[0], soc_converged=False
    )
    unconverged_dossier = replace(
        base_dossier, receipts=tuple(unconverged_receipts)
    )
    write_molecular_soc_dossier_v231(
        unconverged_directory,
        unconverged_dossier,
        dataset=even_first["dataset"],
        identity=even_first["provider"].molecular_soc_contract.identity,
        overwrite=True,
    )
    unconverged_report = audit_molecular_soc_provider_v231(
        even_first["provider"],
        even_first["center"],
        unconverged_directory,
        requirement="protocol",
    )

    tampered_directory = work_directory / "tampered_evidence"
    shutil.copytree(even_first["directory"], tampered_directory)
    reference = base_dossier.evidence.reference
    tampered_reference = replace(
        reference,
        computed_values=(complex(reference.computed_values[0]) + 3.0e-6,),
    )
    tampered_evidence = replace(base_dossier.evidence, reference=tampered_reference)
    tampered_dossier = replace(base_dossier, evidence=tampered_evidence)
    write_molecular_soc_dossier_v231(
        tampered_directory,
        tampered_dossier,
        dataset=even_first["dataset"],
        identity=even_first["provider"].molecular_soc_contract.identity,
        overwrite=True,
    )
    tampered_report = audit_molecular_soc_provider_v231(
        even_first["provider"],
        even_first["center"],
        tampered_directory,
        requirement="protocol",
    )

    disconnected_tracking = replace(
        base_dossier.evidence.tracking,
        record_edges=((0, 1), (1, 2), (2, 3), (5, 6), (6, 7), (7, 8)),
    )
    disconnected_tracking_rejected = _expect_value_error_v231(
        lambda: disconnected_tracking.validate(nrecord=9, nstate=4), "connected"
    )
    traversal = replace(
        base_dossier.artifacts[0], relative_path="../outside.json"
    )
    traversal_rejected = _expect_value_error_v231(traversal.validate, "inside")
    duplicate_receipts = list(base_dossier.receipts)
    duplicate_receipts[1] = replace(
        duplicate_receipts[1],
        output_artifact=duplicate_receipts[0].output_artifact,
    )
    duplicate_output_rejected = _expect_value_error_v231(
        lambda: replace(base_dossier, receipts=tuple(duplicate_receipts)).validate(
            bundle_directory=even_first["directory"]
        ),
        "distinct",
    )
    wrong_environment = replace(
        even_first["provider"].molecular_soc_contract.identity,
        environment_sha256="0" * 64,
    )
    environment_mismatch_rejected = _expect_value_error_v231(
        lambda: base_dossier.validate(
            bundle_directory=even_first["directory"], identity=wrong_environment
        ),
        "environment",
    )
    coordinate_receipts = list(base_dossier.receipts)
    coordinate_receipts[0] = replace(
        coordinate_receipts[0],
        q_bohr=(coordinate_receipts[0].q_bohr[0] + 1.0e-3,),
    )
    coordinate_mismatch_rejected = _expect_value_error_v231(
        lambda: replace(base_dossier, receipts=tuple(coordinate_receipts)).validate(
            bundle_directory=even_first["directory"],
            dataset=even_first["dataset"],
        ),
        "coordinate",
    )

    external = build_v231_admission_bundle(
        work_directory / "synthetic_external_relabel",
        source_kind="external_ab_initio_snapshot",
    )
    inherited_external = audit_molecular_soc_provider_v230(
        external["provider"], external["center"], require_real_backend=True
    )
    external_report = audit_molecular_soc_provider_v231(
        external["provider"],
        external["center"],
        external["dossier_path"],
        requirement="external",
    )
    live = build_v231_admission_bundle(
        work_directory / "synthetic_live_relabel", source_kind="live_ab_initio"
    )
    live_report = audit_molecular_soc_provider_v231(
        live["provider"], live["center"], live["dossier_path"], requirement="live"
    )
    pyscf_probe = probe_pyscf_soc_adapter_v231()
    pyscf_fail_closed = False
    try:
        require_pyscf_soc_adapter_v231()
    except (ImportError, RuntimeError):
        pyscf_fail_closed = True

    return {
        "even_report": even_report.as_dict(),
        "odd_report": odd_report.as_dict(),
        "even_real_report": even_real.as_dict(),
        "odd_real_report": odd_real.as_dict(),
        "dossier_bytes_identical": (
            even_first["dossier_path"].read_bytes()
            == even_second["dossier_path"].read_bytes()
        ),
        "dossier_fingerprints_identical": (
            even_first["dossier"].fingerprint()
            == even_second["dossier"].fingerprint()
        ),
        "raw_artifact_inventories_identical": first_inventory == second_inventory,
        "artifact_corruption_rejected": artifact_corruption_rejected,
        "dossier_corruption_rejected": dossier_corruption_rejected,
        "missing_receipt_rejected": missing_receipt_rejected,
        "unconverged_report": unconverged_report.as_dict(),
        "tampered_evidence_report": tampered_report.as_dict(),
        "disconnected_tracking_rejected": disconnected_tracking_rejected,
        "path_traversal_rejected": traversal_rejected,
        "duplicate_output_rejected": duplicate_output_rejected,
        "environment_mismatch_rejected": environment_mismatch_rejected,
        "coordinate_mismatch_rejected": coordinate_mismatch_rejected,
        "inherited_external_report": inherited_external.as_dict(),
        "external_report": external_report.as_dict(),
        "live_report": live_report.as_dict(),
        "pyscf": pyscf_probe.as_dict(),
        "pyscf_fail_closed": pyscf_fail_closed,
        "even_bundle": even_first,
        "odd_bundle": odd,
    }


def run_v0231_release_benchmark():
    inherited = run_v0230_release_benchmark()
    thresholds = V231AcceptanceThresholds()
    with tempfile.TemporaryDirectory(prefix="gnd-v231-") as temporary:
        controls = _campaign_controls_v231(temporary)
        even = controls["even_report"]
        odd = controls["odd_report"]
        even_real = controls["even_real_report"]
        odd_real = controls["odd_real_report"]
        dossier_checks = even["dossier_protocol_checks"]
        new_checks = {
            "dossier::canonical_bytes_deterministic": controls[
                "dossier_bytes_identical"
            ],
            "dossier::fingerprint_deterministic": controls[
                "dossier_fingerprints_identical"
            ],
            "dossier::raw_artifact_inventory_deterministic": controls[
                "raw_artifact_inventories_identical"
            ],
            "protocol::even_singlet_triplet_passes": even["protocol_passed"],
            "protocol::odd_doublet_passes": odd["protocol_passed"],
            "claims::even_fixture_not_real": not even_real["real_backend_admitted"],
            "claims::odd_fixture_not_real": not odd_real["real_backend_admitted"],
            "binding::replay_dataset": dossier_checks["replay_dataset_binding"],
            "receipts::trajectory_coverage": dossier_checks[
                "trajectory_receipt_coverage"
            ],
            "receipts::all_stages_converged": dossier_checks[
                "all_evidence_calculations_converged"
            ],
            "evidence::independent_reference_derived": dossier_checks[
                "independent_reference_derived"
            ],
            "evidence::basis_convergence_derived": dossier_checks[
                "basis_convergence_derived"
            ],
            "evidence::method_convergence_derived": dossier_checks[
                "method_convergence_derived"
            ],
            "evidence::frame_invariance_derived": dossier_checks[
                "translation_rotation_invariance_derived"
            ],
            "evidence::connected_subspace_tracking_derived": dossier_checks[
                "connected_subspace_tracking_derived"
            ],
            "evidence::summary_matches_raw_observations": dossier_checks[
                "derived_evidence_matches_provider_contract"
            ],
            "integrity::raw_artifact_corruption_rejected": controls[
                "artifact_corruption_rejected"
            ],
            "integrity::dossier_corruption_rejected": controls[
                "dossier_corruption_rejected"
            ],
            "receipts::missing_record_rejected": controls[
                "missing_receipt_rejected"
            ],
            "receipts::unconverged_record_rejected": not controls[
                "unconverged_report"
            ]["protocol_passed"],
            "evidence::tampered_summary_rejected": not controls[
                "tampered_evidence_report"
            ]["protocol_passed"],
            "tracking::disconnected_graph_rejected": controls[
                "disconnected_tracking_rejected"
            ],
            "security::artifact_path_traversal_rejected": controls[
                "path_traversal_rejected"
            ],
            "integrity::duplicate_output_rejected": controls[
                "duplicate_output_rejected"
            ],
            "identity::environment_mismatch_rejected": controls[
                "environment_mismatch_rejected"
            ],
            "identity::receipt_coordinate_mismatch_rejected": controls[
                "coordinate_mismatch_rejected"
            ],
            "admission::v230_summary_alone_is_insufficient": (
                controls["inherited_external_report"]["real_backend_admitted"]
                and not controls["external_report"]["external_snapshot_admitted"]
            ),
            "admission::external_snapshot_requires_parser_attestation": not controls[
                "external_report"
            ]["external_admission_checks"]["backend_artifact_parser_validated"]
            and not controls["external_report"]["external_admission_checks"][
                "executable_artifact_validation"
            ],
            "admission::live_backend_requires_fresh_runtime": not controls[
                "live_report"
            ]["live_admission_checks"]["fresh_runtime_execution_validated"],
            "pyscf::unavailable_or_incomplete_fails_closed": controls[
                "pyscf_fail_closed"
            ],
        }
    new_checks = {name: bool(value) for name, value in new_checks.items()}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.23.1 campaign must define exactly 30 new gates.")
    inherited_checks = {
        f"inherited_v0230::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.23.1 must inherit exactly 93 v0.23.0 gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.23.1 campaign must define exactly 123 total gates.")
    return {
        "release": "v0.23.1",
        "theme": "raw-evidence admission dossiers and executable backend attestation",
        "dossier_protocol": controls["even_report"],
        "odd_doublet_dossier_protocol": controls["odd_report"],
        "negative_controls": {
            key: value
            for key, value in controls.items()
            if key not in {"even_bundle", "odd_bundle"}
        },
        "pyscf": controls["pyscf"],
        "claims": {
            "raw_evidence_admission_protocol_validated": True,
            "calculation_receipt_integrity_validated": True,
            "external_molecular_SOC_snapshot_admitted": False,
            "live_molecular_SOC_backend_admitted": False,
            "ab_initio_SOC_validated": False,
            "live_PySCF_SOC_runtime_validated": False,
            "physical_analytic_SOC_inherited": True,
        },
        "inherited_v0230": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
