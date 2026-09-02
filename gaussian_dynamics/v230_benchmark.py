"""Molecular-SOC admission and deterministic replay campaign for v0.23.0."""

from dataclasses import asdict, dataclass, replace
import importlib.util
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
from .molecular_soc_contract_v230 import (
    MolecularSOCAdmissionContractV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCCapabilitiesV230,
    MolecularSOCValidationEvidenceV230,
    require_trajectory_ready_molecular_soc_v230,
)
from .molecular_soc_replay_v230 import (
    REPLAY_ARRAYS_NAME_V230,
    REPLAY_MANIFEST_NAME_V230,
    FileBackedMolecularSOCProviderV230,
    _canonical_json_bytes_v230,
    _sha256_file_v230,
    _write_deterministic_npz_v230,
    capture_molecular_soc_replay_v230,
    load_molecular_soc_replay_v230,
)
from .pyscf_soc_bridge_v230 import (
    probe_pyscf_soc_runtime_v230,
    require_pyscf_soc_runtime_v230,
)
from .v221_benchmark import run_v0221_release_benchmark


@dataclass(frozen=True)
class V230AcceptanceThresholds:
    roundtrip_matrix_tolerance: float = 0.0
    overlap_roundtrip_tolerance: float = 0.0
    expected_inherited_gates: int = 67
    expected_new_gates: int = 26


def _trajectory_capabilities_v230():
    return MolecularSOCCapabilitiesV230(
        static_soc=True,
        spin_free_derivatives=True,
        soc_derivatives=True,
        derivative_connections=True,
        cross_geometry_overlaps=True,
        deterministic_replay=True,
        analytic_soc_derivatives=True,
    ).validate()


def _fixture_contract_v230(
    *,
    source_kind="validation_fixture",
    all_converged=True,
    evidence=None,
    capabilities=None,
    electron_count=2,
    fixture_name="singlet-triplet",
):
    real_source = source_kind in {
        "external_ab_initio_snapshot",
        "live_ab_initio",
    }
    return MolecularSOCAdmissionContractV230(
        capabilities=(
            _trajectory_capabilities_v230()
            if capabilities is None
            else capabilities
        ),
        identity=MolecularSOCBackendIdentityV230(
            backend_name=(
                f"deterministic analytic {fixture_name} molecular-format fixture"
            ),
            backend_version="1",
            source_kind=source_kind,
            electronic_method=f"analytic {fixture_name} reference",
            basis="closed-form model basis",
            charge=0,
            electron_count=electron_count,
            soc_operator="analytic time-reversal-even SOC",
            scalar_relativistic_method="none",
            derivative_method="analytic component derivatives",
            active_space="complete four-state validation space",
            molecule_name=(
                "negative-control pseudo-molecule"
                if real_source
                else "analytic validation fixture"
            ),
            atom_symbols=(("H", "H") if real_source else ()),
            isotope_masses_amu=((1.007825, 1.007825) if real_source else ()),
            reference_geometry_bohr=(
                ((0.0, 0.0, -0.7), (0.0, 0.0, 0.7)) if real_source else ()
            ),
            calculation_input_sha256=("1" * 64 if real_source else None),
            environment_sha256=("2" * 64 if real_source else None),
            extra={"molecular_accuracy_claim": False},
        ),
        evidence=(
            MolecularSOCValidationEvidenceV230()
            if evidence is None
            else evidence
        ),
        state_tracking_policy="fixed-frame exact overlap fixture",
        coordinate_definition="one explicit generalized coordinate in bohr",
        all_electronic_calculations_converged=all_converged,
    ).validate()


def v230_reference_coordinates(center=0.17):
    center = float(center)
    steps = (1.0e-4, 5.0e-5, 2.5e-5, 1.0e-5)
    values = {center}
    for step in steps:
        values.update((center - step, center + step))
    return np.asarray([[value] for value in sorted(values)], dtype=float)


def build_v230_reference_replay(directory, *, overwrite=False):
    source = AnalyticSingletTripletSOCProviderV220()
    contract = _fixture_contract_v230()
    return capture_molecular_soc_replay_v230(
        directory,
        source,
        v230_reference_coordinates(),
        contract,
        overwrite=overwrite,
    )


def build_v230_doublet_reference_replay(directory, *, overwrite=False):
    source = AnalyticDoubletSOCProviderV220()
    contract = _fixture_contract_v230(
        electron_count=3,
        fixture_name="Kramers-doublet",
    )
    return capture_molecular_soc_replay_v230(
        directory,
        source,
        v230_reference_coordinates(-0.11),
        contract,
        overwrite=overwrite,
    )


def _capture_provider_v230(directory, contract, convergence_flags=None):
    source = AnalyticSingletTripletSOCProviderV220()
    capture_molecular_soc_replay_v230(
        directory,
        source,
        v230_reference_coordinates(),
        contract,
        convergence_flags=convergence_flags,
    )
    return FileBackedMolecularSOCProviderV230(directory)


def _recompute_manifest_fingerprint_v230(path):
    path = Path(path)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    payload = dict(manifest)
    payload.pop("dataset_fingerprint", None)
    manifest["dataset_fingerprint"] = __import__("hashlib").sha256(
        _canonical_json_bytes_v230(payload)
    ).hexdigest()
    path.write_bytes(_canonical_json_bytes_v230(manifest))


def _replay_hardening_v230(work_directory):
    work_directory = Path(work_directory)
    first = work_directory / "first"
    second = work_directory / "second"
    source = AnalyticSingletTripletSOCProviderV220()
    dataset_first = build_v230_reference_replay(first)
    dataset_second = build_v230_reference_replay(second)
    provider = FileBackedMolecularSOCProviderV230(first)
    component_errors = []
    operator_errors = []
    overlap_errors = []
    source_snapshots = []
    replay_snapshots = []
    for q in v230_reference_coordinates():
        source_component = source.components(q)
        replay_component = provider.components(q)
        component_errors.extend(
            float(np.max(np.abs(left - right)))
            for left, right in (
                (source_component.H_spin_free, replay_component.H_spin_free),
                (source_component.K_spin_free, replay_component.K_spin_free),
                (source_component.H_soc, replay_component.H_soc),
                (source_component.K_soc, replay_component.K_soc),
            )
        )
        source_snapshot = source.evaluate_snapshot(q)
        replay_snapshot = provider.evaluate_snapshot(q)
        source_snapshots.append(source_snapshot)
        replay_snapshots.append(replay_snapshot)
        operator_errors.extend(
            float(np.max(np.abs(left - right)))
            for left, right in (
                (source_snapshot.point.H, replay_snapshot.point.H),
                (source_snapshot.point.dH_dq, replay_snapshot.point.dH_dq),
                (source_snapshot.point.connection_q, replay_snapshot.point.connection_q),
                (
                    source_snapshot.point.mass_matrix_q_au,
                    replay_snapshot.point.mass_matrix_q_au,
                ),
            )
        )
    for left in range(len(source_snapshots)):
        for right in range(len(source_snapshots)):
            overlap_errors.append(
                float(
                    np.max(
                        np.abs(
                            source.snapshot_overlap(
                                source_snapshots[left], source_snapshots[right]
                            )
                            - provider.snapshot_overlap(
                                replay_snapshots[left], replay_snapshots[right]
                            )
                        )
                    )
                )
            )

    coordinate_miss_rejected = False
    try:
        provider.evaluate_snapshot(np.asarray([0.17123]))
    except KeyError:
        coordinate_miss_rejected = True

    array_corruption = work_directory / "array_corruption"
    shutil.copytree(first, array_corruption)
    arrays_path = array_corruption / REPLAY_ARRAYS_NAME_V230
    raw = bytearray(arrays_path.read_bytes())
    raw[len(raw) // 2] ^= 0x01
    arrays_path.write_bytes(raw)
    array_corruption_rejected = False
    try:
        load_molecular_soc_replay_v230(array_corruption)
    except ValueError as exc:
        array_corruption_rejected = "integrity" in str(exc)

    manifest_corruption = work_directory / "manifest_corruption"
    shutil.copytree(first, manifest_corruption)
    manifest_path = manifest_corruption / REPLAY_MANIFEST_NAME_V230
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["provenance"]["model_name"] = "changed without fingerprint update"
    manifest_path.write_bytes(_canonical_json_bytes_v230(manifest))
    manifest_corruption_rejected = False
    try:
        load_molecular_soc_replay_v230(manifest_corruption)
    except ValueError as exc:
        manifest_corruption_rejected = "manifest integrity" in str(exc)

    overlap_corruption = work_directory / "overlap_corruption"
    shutil.copytree(first, overlap_corruption)
    overlap_arrays_path = overlap_corruption / REPLAY_ARRAYS_NAME_V230
    with np.load(overlap_arrays_path, allow_pickle=False) as archive:
        arrays = {name: archive[name].copy() for name in archive.files}
    arrays["overlaps"][0, 1, 0, 0] += 0.1
    _write_deterministic_npz_v230(overlap_arrays_path, arrays)
    overlap_manifest_path = overlap_corruption / REPLAY_MANIFEST_NAME_V230
    overlap_manifest = json.loads(overlap_manifest_path.read_text(encoding="utf-8"))
    overlap_manifest["arrays_sha256"] = _sha256_file_v230(overlap_arrays_path)
    overlap_manifest_path.write_bytes(_canonical_json_bytes_v230(overlap_manifest))
    _recompute_manifest_fingerprint_v230(overlap_manifest_path)
    overlap_corruption_rejected = False
    try:
        load_molecular_soc_replay_v230(overlap_corruption)
    except ValueError as exc:
        overlap_corruption_rejected = "overlaps violate" in str(exc)

    return {
        "dataset_fingerprint": dataset_first.dataset_fingerprint,
        "record_count": int(len(dataset_first.q)),
        "maximum_component_roundtrip_error": max(component_errors, default=0.0),
        "maximum_operator_roundtrip_error": max(operator_errors, default=0.0),
        "maximum_overlap_roundtrip_error": max(overlap_errors, default=0.0),
        "manifest_bytes_identical": (
            (first / REPLAY_MANIFEST_NAME_V230).read_bytes()
            == (second / REPLAY_MANIFEST_NAME_V230).read_bytes()
        ),
        "array_bytes_identical": (
            (first / REPLAY_ARRAYS_NAME_V230).read_bytes()
            == (second / REPLAY_ARRAYS_NAME_V230).read_bytes()
        ),
        "dataset_fingerprints_identical": (
            dataset_first.dataset_fingerprint == dataset_second.dataset_fingerprint
        ),
        "coordinate_miss_rejected": coordinate_miss_rejected,
        "array_corruption_rejected": array_corruption_rejected,
        "manifest_corruption_rejected": manifest_corruption_rejected,
        "overlap_corruption_rejected": overlap_corruption_rejected,
        "provider": provider,
    }


def _odd_replay_v230(work_directory):
    work_directory = Path(work_directory)
    source = AnalyticDoubletSOCProviderV220()
    dataset = build_v230_doublet_reference_replay(work_directory)
    provider = FileBackedMolecularSOCProviderV230(work_directory)
    component_errors = []
    overlap_errors = []
    source_snapshots = []
    replay_snapshots = []
    for q in v230_reference_coordinates(-0.11):
        expected = source.components(q)
        actual = provider.components(q)
        component_errors.extend(
            float(np.max(np.abs(left - right)))
            for left, right in (
                (expected.H_spin_free, actual.H_spin_free),
                (expected.K_spin_free, actual.K_spin_free),
                (expected.H_soc, actual.H_soc),
                (expected.K_soc, actual.K_soc),
            )
        )
        source_snapshots.append(source.evaluate_snapshot(q))
        replay_snapshots.append(provider.evaluate_snapshot(q))
    for left in range(len(source_snapshots)):
        for right in range(len(source_snapshots)):
            overlap_errors.append(
                float(
                    np.max(
                        np.abs(
                            source.snapshot_overlap(
                                source_snapshots[left], source_snapshots[right]
                            )
                            - provider.snapshot_overlap(
                                replay_snapshots[left], replay_snapshots[right]
                            )
                        )
                    )
                )
            )
    protocol = audit_molecular_soc_provider_v230(
        provider, np.asarray([-0.11]), require_real_backend=False
    )
    real = audit_molecular_soc_provider_v230(
        provider, np.asarray([-0.11]), require_real_backend=True
    )
    return {
        "dataset_fingerprint": dataset.dataset_fingerprint,
        "record_count": int(len(dataset.q)),
        "electron_parity": provider.soc_symmetry_contract.electron_parity,
        "maximum_component_roundtrip_error": max(component_errors, default=0.0),
        "maximum_overlap_roundtrip_error": max(overlap_errors, default=0.0),
        "protocol_report": protocol.as_dict(),
        "real_admission_report": real.as_dict(),
    }


class _StaticContractProbeV230:
    def __init__(self, symmetry):
        self.soc_symmetry_contract = symmetry
        self.molecular_soc_contract = _fixture_contract_v230(
            capabilities=MolecularSOCCapabilitiesV230(
                static_soc=True,
                deterministic_replay=True,
            )
        )


def _admission_hardening_v230(work_directory, reference_provider):
    work_directory = Path(work_directory)
    center = np.asarray([0.17])
    protocol = audit_molecular_soc_provider_v230(
        reference_provider, center, require_real_backend=False
    )
    real = audit_molecular_soc_provider_v230(
        reference_provider, center, require_real_backend=True
    )

    static_probe = _StaticContractProbeV230(reference_provider.soc_symmetry_contract)
    static_rejected = False
    try:
        require_trajectory_ready_molecular_soc_v230(static_probe)
    except ValueError:
        static_rejected = True

    unconverged_flags = np.ones(len(v230_reference_coordinates()), dtype=bool)
    unconverged_flags[-1] = False
    unconverged_contract = _fixture_contract_v230(all_converged=False)
    unconverged_provider = _capture_provider_v230(
        work_directory / "unconverged",
        unconverged_contract,
        convergence_flags=unconverged_flags,
    )
    unconverged = audit_molecular_soc_provider_v230(
        unconverged_provider, center, require_real_backend=False
    )

    reference = dict(
        independent_reference_id="negative-control reference",
        independent_reference_error=1.0e-6,
        independent_reference_tolerance=1.0e-5,
    )
    valid_basis = dict(
        basis_levels=("small", "medium", "large"),
        basis_changes=(2.0e-4, 5.0e-6),
        basis_tolerance=1.0e-5,
    )
    valid_method = dict(
        method_levels=("method-a", "method-b"),
        method_changes=(4.0e-6,),
        method_tolerance=1.0e-5,
    )
    valid_invariance = dict(
        translation_residual=1.0e-9,
        rotation_residual=2.0e-9,
        frame_invariance_tolerance=1.0e-8,
    )
    valid_tracking = dict(
        tracking_minimum_overlap=0.92,
        tracking_minimum_margin=0.18,
        tracking_overlap_threshold=0.80,
        tracking_margin_threshold=0.10,
    )
    evidence_cases = {
        "missing_reference": MolecularSOCValidationEvidenceV230(
            **valid_basis, **valid_method, **valid_invariance, **valid_tracking
        ),
        "missing_basis": MolecularSOCValidationEvidenceV230(
            **reference, **valid_method, **valid_invariance, **valid_tracking
        ),
        "missing_method": MolecularSOCValidationEvidenceV230(
            **reference, **valid_basis, **valid_invariance, **valid_tracking
        ),
        "missing_invariance": MolecularSOCValidationEvidenceV230(
            **reference, **valid_basis, **valid_method, **valid_tracking
        ),
        "missing_tracking": MolecularSOCValidationEvidenceV230(
            **reference, **valid_basis, **valid_method, **valid_invariance
        ),
    }
    evidence_reports = {}
    for name, evidence in evidence_cases.items():
        contract = _fixture_contract_v230(
            source_kind="external_ab_initio_snapshot",
            evidence=evidence,
        )
        provider = _capture_provider_v230(
            work_directory / name,
            contract,
            convergence_flags=np.ones(len(v230_reference_coordinates()), dtype=bool),
        )
        evidence_reports[name] = audit_molecular_soc_provider_v230(
            provider, center, require_real_backend=True
        ).as_dict()

    return {
        "protocol_report": protocol.as_dict(),
        "real_admission_report": real.as_dict(),
        "static_tier": static_probe.molecular_soc_contract.capabilities.as_dict(),
        "static_rejected_for_nuclear_motion": static_rejected,
        "unconverged_report": unconverged.as_dict(),
        "evidence_negative_controls": evidence_reports,
    }


def run_v0230_release_benchmark():
    inherited = run_v0221_release_benchmark()
    thresholds = V230AcceptanceThresholds()
    with tempfile.TemporaryDirectory(prefix="gnd-v230-") as temporary:
        replay = _replay_hardening_v230(Path(temporary) / "replay")
        provider = replay.pop("provider")
        admission = _admission_hardening_v230(Path(temporary) / "admission", provider)
        odd_replay = _odd_replay_v230(Path(temporary) / "odd_replay")
    pyscf_probe = probe_pyscf_soc_runtime_v230()
    pyscf_fail_closed = False
    if not pyscf_probe.installed:
        try:
            require_pyscf_soc_runtime_v230()
        except ImportError:
            pyscf_fail_closed = True
    else:
        pyscf_fail_closed = not pyscf_probe.live_soc_adapter_validated

    protocol = admission["protocol_report"]
    real = admission["real_admission_report"]
    unconverged = admission["unconverged_report"]
    evidence = admission["evidence_negative_controls"]
    new_checks = {
        "capability::static_tier_declared": (
            admission["static_tier"]["tier"] == "static_soc"
        ),
        "capability::static_rejected_for_nuclear_motion": admission[
            "static_rejected_for_nuclear_motion"
        ],
        "capability::trajectory_tier_complete": (
            protocol["capability_tier"] == "trajectory_ready"
        ),
        "replay::component_roundtrip_exact": (
            replay["maximum_component_roundtrip_error"]
            <= thresholds.roundtrip_matrix_tolerance
        ),
        "replay::operator_roundtrip_exact": (
            replay["maximum_operator_roundtrip_error"]
            <= thresholds.roundtrip_matrix_tolerance
        ),
        "replay::overlap_roundtrip_exact": (
            replay["maximum_overlap_roundtrip_error"]
            <= thresholds.overlap_roundtrip_tolerance
        ),
        "replay::manifest_deterministic": replay["manifest_bytes_identical"],
        "replay::arrays_deterministic": replay["array_bytes_identical"],
        "replay::dataset_fingerprint_deterministic": replay[
            "dataset_fingerprints_identical"
        ],
        "replay::coordinate_miss_rejected": replay["coordinate_miss_rejected"],
        "replay::array_corruption_rejected": replay["array_corruption_rejected"],
        "replay::manifest_corruption_rejected": replay[
            "manifest_corruption_rejected"
        ],
        "replay::overlap_corruption_rejected": replay[
            "overlap_corruption_rejected"
        ],
        "admission::fixture_protocol_passes": protocol["protocol_passed"],
        "admission::fixture_not_real_backend": (
            not real["real_backend_admitted"]
            and not real["real_admission_checks"]["real_ab_initio_source"]
        ),
        "admission::unconverged_records_rejected": (
            not unconverged["protocol_passed"]
            and not unconverged["protocol_checks"][
                "all_electronic_calculations_converged"
            ]
        ),
        "admission::missing_reference_rejected": (
            not evidence["missing_reference"]["real_backend_admitted"]
            and not evidence["missing_reference"]["real_admission_checks"][
                "independent_reference_evidence"
            ]
        ),
        "admission::missing_basis_convergence_rejected": (
            not evidence["missing_basis"]["real_backend_admitted"]
            and not evidence["missing_basis"]["real_admission_checks"][
                "basis_convergence_evidence"
            ]
        ),
        "admission::missing_method_convergence_rejected": (
            not evidence["missing_method"]["real_backend_admitted"]
            and not evidence["missing_method"]["real_admission_checks"][
                "method_convergence_evidence"
            ]
        ),
        "admission::missing_frame_invariance_rejected": (
            not evidence["missing_invariance"]["real_backend_admitted"]
            and not evidence["missing_invariance"]["real_admission_checks"][
                "translation_rotation_invariance"
            ]
        ),
        "admission::missing_tracking_quality_rejected": (
            not evidence["missing_tracking"]["real_backend_admitted"]
            and not evidence["missing_tracking"]["real_admission_checks"][
                "state_tracking_quality"
            ]
        ),
        "admission::symmetry_and_component_derivatives": (
            protocol["protocol_checks"]["single_electron_parity_and_charge"]
            and protocol["protocol_checks"]["component_resolved_derivatives"]
            and protocol["protocol_checks"]["cross_geometry_differentials"]
        ),
        "replay::odd_doublet_components_exact": (
            odd_replay["electron_parity"] == "odd"
            and odd_replay["maximum_component_roundtrip_error"]
            <= thresholds.roundtrip_matrix_tolerance
            and odd_replay["maximum_overlap_roundtrip_error"]
            <= thresholds.overlap_roundtrip_tolerance
        ),
        "admission::odd_doublet_protocol_passes": odd_replay[
            "protocol_report"
        ]["protocol_passed"],
        "admission::odd_doublet_fixture_not_real": (
            not odd_replay["real_admission_report"]["real_backend_admitted"]
            and not odd_replay["real_admission_report"]["real_admission_checks"][
                "real_ab_initio_source"
            ]
        ),
        "pyscf::unavailable_or_unvalidated_fails_closed": pyscf_fail_closed,
    }
    new_checks = {name: bool(value) for name, value in new_checks.items()}
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError("v0.23.0 campaign must define exactly 26 new gates.")
    inherited_checks = {
        f"inherited_v0221::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.23.0 must inherit exactly 67 v0.22.1 gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != 93:
        raise AssertionError("v0.23.0 campaign must define exactly 93 total gates.")
    return {
        "release": "v0.23.0",
        "theme": "molecular SOC backend admission protocol and deterministic replay",
        "replay": replay,
        "admission": admission,
        "odd_doublet_replay": odd_replay,
        "pyscf": pyscf_probe.as_dict(),
        "claims": {
            "molecular_SOC_protocol_validated": True,
            "deterministic_replay_validated": True,
            "real_molecular_SOC_backend_admitted": False,
            "ab_initio_SOC_validated": False,
            "live_PySCF_SOC_runtime_validated": False,
            "physical_analytic_SOC_inherited": True,
        },
        "inherited_v0221": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
