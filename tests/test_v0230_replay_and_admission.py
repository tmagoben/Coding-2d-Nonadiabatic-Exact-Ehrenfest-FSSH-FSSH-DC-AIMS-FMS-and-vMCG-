from pathlib import Path

import numpy as np
import pytest

from gaussian_dynamics import (
    FileBackedMolecularSOCProviderV230,
    audit_molecular_soc_provider_v230,
    load_molecular_soc_replay_v230,
    probe_pyscf_soc_runtime_v230,
    require_pyscf_soc_runtime_v230,
)
from gaussian_dynamics.v230_benchmark import (
    _admission_hardening_v230,
    _replay_hardening_v230,
    build_v230_reference_replay,
    v230_reference_coordinates,
)


def test_reference_replay_roundtrip_is_exact_and_deterministic(tmp_path):
    result = _replay_hardening_v230(tmp_path)

    assert result["record_count"] == 9
    assert result["maximum_component_roundtrip_error"] == 0.0
    assert result["maximum_operator_roundtrip_error"] == 0.0
    assert result["maximum_overlap_roundtrip_error"] == 0.0
    assert result["manifest_bytes_identical"]
    assert result["array_bytes_identical"]
    assert result["dataset_fingerprints_identical"]


def test_replay_rejects_unknown_coordinate_without_interpolation(tmp_path):
    build_v230_reference_replay(tmp_path)
    provider = FileBackedMolecularSOCProviderV230(tmp_path)

    with pytest.raises(KeyError, match="interpolation and extrapolation are forbidden"):
        provider.evaluate_snapshot(np.asarray([0.17123]))


def test_replay_detects_array_manifest_and_overlap_corruption(tmp_path):
    result = _replay_hardening_v230(tmp_path)

    assert result["array_corruption_rejected"]
    assert result["manifest_corruption_rejected"]
    assert result["overlap_corruption_rejected"]


def test_loaded_dataset_retains_contract_and_convergence_identity(tmp_path):
    captured = build_v230_reference_replay(tmp_path)
    loaded = load_molecular_soc_replay_v230(tmp_path)

    assert captured.dataset_fingerprint == loaded.dataset_fingerprint
    assert loaded.molecular_soc_contract.capabilities.trajectory_ready
    assert loaded.molecular_soc_contract.identity.source_kind == "validation_fixture"
    assert np.all(loaded.converged)


def test_fixture_passes_protocol_but_cannot_be_admitted_as_real(tmp_path):
    build_v230_reference_replay(tmp_path)
    provider = FileBackedMolecularSOCProviderV230(tmp_path)
    protocol = audit_molecular_soc_provider_v230(
        provider, np.asarray([0.17]), require_real_backend=False
    )
    real = audit_molecular_soc_provider_v230(
        provider, np.asarray([0.17]), require_real_backend=True
    )

    assert protocol.protocol_passed
    assert protocol.passed
    assert not real.real_backend_admitted
    assert not real.passed
    assert not real.real_admission_checks["real_ab_initio_source"]


def test_unconverged_and_missing_evidence_controls_are_independent(tmp_path):
    build_v230_reference_replay(tmp_path / "reference")
    reference = FileBackedMolecularSOCProviderV230(tmp_path / "reference")
    result = _admission_hardening_v230(tmp_path / "controls", reference)

    assert not result["unconverged_report"]["protocol_passed"]
    controls = result["evidence_negative_controls"]
    assert not controls["missing_reference"]["real_admission_checks"][
        "independent_reference_evidence"
    ]
    assert not controls["missing_basis"]["real_admission_checks"][
        "basis_convergence_evidence"
    ]
    assert not controls["missing_method"]["real_admission_checks"][
        "method_convergence_evidence"
    ]
    assert not controls["missing_invariance"]["real_admission_checks"][
        "translation_rotation_invariance"
    ]
    assert not controls["missing_tracking"]["real_admission_checks"][
        "state_tracking_quality"
    ]


def test_replay_snapshot_tokens_reject_foreign_dataset(tmp_path):
    build_v230_reference_replay(tmp_path / "first")
    first = FileBackedMolecularSOCProviderV230(tmp_path / "first")
    snapshot = first.evaluate_snapshot(v230_reference_coordinates()[0])
    foreign = FileBackedMolecularSOCProviderV230(tmp_path / "first")
    token = snapshot.wavefunction_snapshot
    object.__setattr__(token, "dataset_fingerprint", "0" * 64)

    with pytest.raises(ValueError, match="cross-dataset"):
        foreign.snapshot_overlap(snapshot, snapshot)


def test_pyscf_runtime_probe_and_requirement_fail_closed_when_unavailable():
    probe = probe_pyscf_soc_runtime_v230()

    if probe.installed:
        assert not probe.live_soc_adapter_validated
    else:
        with pytest.raises(ImportError, match="no live backend"):
            require_pyscf_soc_runtime_v230()
