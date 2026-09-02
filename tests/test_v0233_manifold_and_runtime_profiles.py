import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from gaussian_dynamics.analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
)
from gaussian_dynamics.manifold_transport_v233 import (
    audit_complete_manifold_transport_v233,
)
from gaussian_dynamics.molecular_soc_convention_v233 import (
    analytic_soc_convention_v233,
    audit_molecular_soc_convention_v233,
    require_exact_molecular_soc_convention_v233,
)
from gaussian_dynamics.molecular_soc_replay_v233 import (
    FileBackedMolecularSOCProviderV233,
    migrate_molecular_soc_replay_v230_to_v233,
)
from gaussian_dynamics.nac_compatibility_v233 import (
    LegacyReplayMigrationAttestationV233,
    analytic_nac_convention_v233,
)
from gaussian_dynamics.provider_numerical_identity_v233 import (
    require_provider_numerical_identity_v233,
)
from gaussian_dynamics.runtime_compatibility_v233 import (
    assess_runtime_compatibility_v233,
    release_locked_runtime_profile_v233,
    scientifically_compatible_runtime_profile_v233,
)
from gaussian_dynamics.soc_admission_v221 import SOCSymmetryContractV221
from gaussian_dynamics.v230_benchmark import build_v230_reference_replay


def _random_unitary(seed, dimension):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
        size=(dimension, dimension)
    )
    q, _ = np.linalg.qr(raw)
    return q


def _block_unitary(seed):
    matrix = np.zeros((4, 4), dtype=complex)
    matrix[:2, :2] = _random_unitary(seed, 2)
    matrix[2:, 2:] = _random_unitary(seed + 1, 2)
    return matrix


def _transform_symmetry(contract, gauge):
    return SOCSymmetryContractV221(
        electron_parity=contract.electron_parity,
        time_reversal_matrix=gauge.conj().T @ contract.time_reversal_matrix @ gauge.conj(),
        projectors={
            name: gauge.conj().T @ projector @ gauge
            for name, projector in contract.projectors.items()
        },
        external_magnetic_field=contract.external_magnetic_field,
    )


def test_arbitrary_endpoint_gauges_preserve_complete_kramers_manifolds():
    provider = AnalyticDoubletSOCProviderV220()
    base = provider.soc_symmetry_contract
    left_gauge = _block_unitary(2331)
    right_gauge = _block_unitary(2333)
    overlap = left_gauge.conj().T @ right_gauge
    report = audit_complete_manifold_transport_v233(
        overlap,
        provider.provenance.model_space,
        _transform_symmetry(base, left_gauge),
        right_symmetry_contract=_transform_symmetry(base, right_gauge),
    )
    assert report.passed
    assert report.electron_parity == "odd"
    assert [item.dimension for item in report.manifold_blocks] == [2, 2]
    assert report.time_reversal_covariance_residual < 1e-10


def test_singlet_triplet_complete_manifolds_pass_and_cross_leakage_fails():
    provider = AnalyticSingletTripletSOCProviderV220()
    identity = np.eye(4, dtype=complex)
    passed = audit_complete_manifold_transport_v233(
        identity,
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
    )
    assert passed.passed
    assert sorted(item.dimension for item in passed.manifold_blocks) == [1, 3]

    swapped = identity.copy()
    swapped[:, [0, 1]] = swapped[:, [1, 0]]
    failed = audit_complete_manifold_transport_v233(
        swapped,
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
    )
    assert not failed.passed
    assert not failed.checks["assigned_manifold_retention"]
    assert not failed.checks["competing_manifold_leakage"]


def test_broken_kramers_overlap_is_detected_by_time_reversal_covariance():
    provider = AnalyticDoubletSOCProviderV220()
    overlap = np.diag([1.0, 0.8, 1.0, 1.0]).astype(complex)
    report = audit_complete_manifold_transport_v233(
        overlap,
        provider.provenance.model_space,
        provider.soc_symmetry_contract,
    )
    assert not report.passed
    assert not report.checks["time_reversal_covariance"]


def test_release_lock_and_scientific_compatibility_are_independent():
    root = Path(__file__).resolve().parents[1]
    evidence = json.loads(
        (root / "results/v0232_pyscf_runtime_evidence.json").read_text(
            encoding="utf-8"
        )
    )["runtime"]
    locked = assess_runtime_compatibility_v233(
        evidence, release_locked_runtime_profile_v233()
    )
    compatible = assess_runtime_compatibility_v233(
        evidence, scientifically_compatible_runtime_profile_v233()
    )
    assert locked.compatible
    assert compatible.compatible

    portable = dict(evidence)
    portable["python_executable_sha256"] = "0" * 64
    portable["platform"] = "Linux-different-kernel"
    assert not assess_runtime_compatibility_v233(
        portable, release_locked_runtime_profile_v233()
    ).compatible
    assert assess_runtime_compatibility_v233(
        portable, scientifically_compatible_runtime_profile_v233()
    ).compatible


def test_soc_matrix_convention_binds_state_order_units_and_derivatives():
    for provider, coordinate in (
        (AnalyticSingletTripletSOCProviderV220(), np.asarray([0.17])),
        (AnalyticDoubletSOCProviderV220(), np.asarray([-0.11])),
    ):
        convention = analytic_soc_convention_v233(provider)
        report = audit_molecular_soc_convention_v233(
            provider.components(coordinate),
            provider.provenance,
            provider.soc_symmetry_contract,
            convention,
        )
        assert report.passed

        wrong_order = replace(
            convention, state_order=tuple(reversed(convention.state_order))
        )
        failed = audit_molecular_soc_convention_v233(
            provider.components(coordinate),
            provider.provenance,
            provider.soc_symmetry_contract,
            wrong_order,
        )
        assert not failed.passed
        assert not failed.checks["exact_state_order"]

        wrong_prefactor = replace(
            convention, prefactor_convention="different hidden prefactor"
        )
        try:
            require_exact_molecular_soc_convention_v233(
                wrong_prefactor, convention
            )
        except ValueError as exc:
            assert "differs from the trusted" in str(exc)
        else:
            raise AssertionError("SOC prefactor mismatch crossed the trust anchor")


def test_replay_provider_numerical_identity_quarantines_legacy_provider(tmp_path):
    legacy = build_v230_reference_replay(tmp_path / "legacy")
    convention = analytic_nac_convention_v233()
    migrated = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        tmp_path / "v233",
        nac_convention=convention,
        migration_attestation=LegacyReplayMigrationAttestationV233(
            legacy_dataset_fingerprint=legacy.dataset_fingerprint,
            nac_disposition="not_pyscf_derived",
            evidence="analytic fixture emits internal d[i,j] directly",
        ),
    )
    provider = FileBackedMolecularSOCProviderV233(migrated.manifest_path)
    identity = provider.numerical_identity_v233
    assert require_provider_numerical_identity_v233(provider, identity) == identity

    class LegacyProvider:
        provenance = provider.provenance

    try:
        require_provider_numerical_identity_v233(LegacyProvider(), identity)
    except ValueError as exc:
        assert "legacy providers are quarantined" in str(exc)
    else:
        raise AssertionError("legacy provider crossed the numerical identity gate")
