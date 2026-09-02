"""Finite-manifold transport and compatibility campaign for v0.23.3."""

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import tempfile

import numpy as np

from .analytic_soc_models_v220 import (
    AnalyticDoubletSOCProviderV220,
    AnalyticSingletTripletSOCProviderV220,
)
from .finite_manifold_transport_v233 import (
    FiniteManifoldOverlapPolicyV233,
    analyze_finite_manifold_overlap_v233,
    certified_transport_from_overlap_v233,
    certify_reciprocal_transport_pair_v233,
)
from .manifold_transport_v233 import audit_complete_manifold_transport_v233
from .molecular_soc_convention_v233 import (
    analytic_soc_convention_v233,
    audit_molecular_soc_convention_v233,
    require_exact_molecular_soc_convention_v233,
)
from .molecular_soc_replay_v233 import (
    FileBackedMolecularSOCProviderV233,
    capture_molecular_soc_replay_v233,
    migrate_molecular_soc_replay_v230_to_v233,
)
from .nac_compatibility_v233 import (
    LegacyReplayMigrationAttestationV233,
    analytic_nac_convention_v233,
    corrected_pyscf_nac_convention_v233,
    require_snapshot_nac_identity_v233,
)
from .overlap_transport import nearest_unitary
from .provider_numerical_identity_v233 import (
    build_provider_numerical_identity_v233,
    require_provider_numerical_identity_v233,
)
from .runtime_compatibility_v233 import (
    assess_runtime_compatibility_v233,
    release_locked_runtime_profile_v233,
    scientifically_compatible_runtime_profile_v233,
)
from .soc_admission_v221 import SOCSymmetryContractV221
from .v230_benchmark import (
    _fixture_contract_v230,
    build_v230_reference_replay,
    v230_reference_coordinates,
)
from .v232_benchmark import run_v0232_release_benchmark


@dataclass(frozen=True)
class V233AcceptanceThresholds:
    expected_inherited_gates: int = 168
    expected_new_gates: int = 40
    expected_total_gates: int = 208


def _raises_v233(callable_object, exceptions, text=None):
    try:
        callable_object()
    except exceptions as exc:
        return text is None or text in str(exc)
    return False


def _random_unitary_v233(seed, dimension):
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(dimension, dimension)) + 1j * rng.normal(
        size=(dimension, dimension)
    )
    q, r = np.linalg.qr(raw)
    diagonal = np.diag(r)
    phases = np.where(np.abs(diagonal) == 0.0, 1.0, diagonal / np.abs(diagonal))
    return q @ np.diag(phases.conj())


def _block_unitary_v233(seed):
    matrix = np.zeros((4, 4), dtype=complex)
    matrix[:2, :2] = _random_unitary_v233(seed, 2)
    matrix[2:, 2:] = _random_unitary_v233(seed + 1, 2)
    return matrix


def _transform_symmetry_v233(contract, gauge):
    return SOCSymmetryContractV221(
        electron_parity=contract.electron_parity,
        time_reversal_matrix=(
            gauge.conj().T @ contract.time_reversal_matrix @ gauge.conj()
        ),
        projectors={
            name: gauge.conj().T @ projector @ gauge
            for name, projector in contract.projectors.items()
        },
        external_magnetic_field=contract.external_magnetic_field,
    )


def _transport_controls_v233():
    contraction = np.diag([0.92, 0.74, 0.61]).astype(complex)
    result = certified_transport_from_overlap_v233(contraction)
    left = _random_unitary_v233(23301, 3)
    right = _random_unitary_v233(23302, 3)
    transformed = certified_transport_from_overlap_v233(
        left.conj().T @ contraction @ right
    )
    expected = left.conj().T @ result.right_to_left_transport @ right
    strict_policy = FiniteManifoldOverlapPolicyV233(
        minimum_retained_singular_value=0.5,
        maximum_condition_number=10.0,
    )
    low_retention = analyze_finite_manifold_overlap_v233(
        np.diag([0.9, 0.01]), policy=strict_policy
    )
    pair = certify_reciprocal_transport_pair_v233(
        left.conj().T @ contraction @ right,
        right.conj().T @ contraction @ left,
    )
    corrupted_reverse = right.conj().T @ contraction @ left
    corrupted_reverse = corrupted_reverse.copy()
    corrupted_reverse[0, 0] += 1.0e-3
    return {
        "raw_contraction_physically_consistent": result.physically_consistent,
        "raw_contraction_is_nonunitary": not np.allclose(
            contraction.conj().T @ contraction, np.eye(3)
        ),
        "polar_transport_is_unitary": np.allclose(
            result.right_to_left_transport.conj().T
            @ result.right_to_left_transport,
            np.eye(3),
            atol=1.0e-12,
        ),
        "raw_overlap_and_transport_are_distinct": not np.allclose(
            result.overlap, result.right_to_left_transport
        ),
        "polar_transport_is_gauge_covariant": np.allclose(
            transformed.right_to_left_transport, expected, atol=1.0e-12
        ),
        "spectral_expansion_rejected": _raises_v233(
            lambda: certified_transport_from_overlap_v233(
                np.diag([1.01, 0.9])
            ),
            (ValueError,),
            "physically inconsistent",
        ),
        "low_retention_fails_trajectory_readiness": bool(
            low_retention.physically_consistent
            and not low_retention.trajectory_ready
        ),
        "ill_conditioning_reported_independently": (
            "ill_conditioned_overlap" in low_retention.failed_quality_checks
        ),
        "raw_and_transport_reciprocity_certified": bool(
            pair.overlap_reciprocity_residual < 1.0e-12
            and pair.transport_reciprocity_residual < 1.0e-12
        ),
        "reciprocity_corruption_rejected": _raises_v233(
            lambda: certify_reciprocal_transport_pair_v233(
                left.conj().T @ contraction @ right, corrupted_reverse
            ),
            (ValueError,),
            "adjoint reciprocity",
        ),
        "rank_loss_rejected_by_transport_consumer": _raises_v233(
            lambda: nearest_unitary(np.diag([1.0, 0.0])),
            (ValueError,),
            "not trajectory ready",
        ),
        "diagnostics": result.as_dict(),
        "low_retention_diagnostics": low_retention.as_dict(),
    }


def _replay_and_nac_controls_v233(directory):
    directory = Path(directory)
    legacy = build_v230_reference_replay(directory / "legacy")
    convention = analytic_nac_convention_v233()
    attestation = LegacyReplayMigrationAttestationV233(
        legacy_dataset_fingerprint=legacy.dataset_fingerprint,
        nac_disposition="not_pyscf_derived",
        evidence="analytic fixture emits internal d[i,j] directly",
    )
    first = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        directory / "first",
        nac_convention=convention,
        migration_attestation=attestation,
    )
    second = migrate_molecular_soc_replay_v230_to_v233(
        legacy.manifest_path,
        directory / "second",
        nac_convention=convention,
        migration_attestation=attestation,
    )
    provider = FileBackedMolecularSOCProviderV233(first.manifest_path)
    coordinates = v230_reference_coordinates()
    left = provider.evaluate_snapshot(coordinates[0])
    right = provider.evaluate_snapshot(coordinates[1])
    raw = provider.snapshot_overlap(left, right)
    transport = provider.snapshot_transport(left, right)
    stored_transport = np.asarray(first.overlap_transports, dtype=complex).copy()
    stored_transport[0, 1] = 0.9 * np.eye(stored_transport.shape[-1])
    raw_substitution_rejected = _raises_v233(
        lambda: replace(first, overlap_transports=stored_transport).validate(),
        (ValueError,),
        "transport differs",
    )
    unknown = replace(attestation, nac_disposition="unknown")
    wrong_sign = replace(attestation, nac_disposition="requires_sign_correction")
    full_pyscf = corrected_pyscf_nac_convention_v233()
    etf_pyscf = corrected_pyscf_nac_convention_v233(use_etfs=True)

    direct_source = AnalyticSingletTripletSOCProviderV220()
    direct = capture_molecular_soc_replay_v233(
        directory / "direct",
        direct_source,
        coordinates,
        _fixture_contract_v230(),
        nac_convention=convention,
    )
    identity = provider.numerical_identity_v233
    changed_nac_identity = build_provider_numerical_identity_v233(
        provider.provenance,
        corrected_pyscf_nac_convention_v233(),
        overlap_policy=provider.dataset.overlap_policy,
    )
    changed_policy = build_provider_numerical_identity_v233(
        provider.provenance,
        convention,
        overlap_policy=FiniteManifoldOverlapPolicyV233(
            minimum_retained_singular_value=0.6
        ),
    )
    class LegacyProvider:
        provenance = provider.provenance

    return {
        "format_version_two": json.loads(
            first.manifest_path.read_text(encoding="utf-8")
        )["format_version"]
        == 2,
        "migration_manifest_deterministic": (
            first.manifest_path.read_bytes() == second.manifest_path.read_bytes()
        ),
        "migration_arrays_deterministic": (
            first.arrays_path.read_bytes() == second.arrays_path.read_bytes()
        ),
        "migration_fingerprint_deterministic": (
            first.dataset_fingerprint == second.dataset_fingerprint
        ),
        "attestation_bound_to_legacy_fingerprint": (
            first.migration_attestation["legacy_dataset_fingerprint"]
            == legacy.dataset_fingerprint
        ),
        "provider_raw_overlap_roundtrip": np.array_equal(
            raw, first.overlaps[0, 1]
        ),
        "provider_transport_roundtrip": np.array_equal(
            transport, first.overlap_transports[0, 1]
        ),
        "raw_overlap_substitution_rejected": raw_substitution_rejected,
        "unknown_legacy_nac_quarantined": _raises_v233(
            lambda: migrate_molecular_soc_replay_v230_to_v233(
                legacy.manifest_path,
                directory / "unknown",
                nac_convention=convention,
                migration_attestation=unknown,
            ),
            (ValueError,),
            "quarantined",
        ),
        "wrong_sign_legacy_nac_quarantined": _raises_v233(
            lambda: migrate_molecular_soc_replay_v230_to_v233(
                legacy.manifest_path,
                directory / "wrong_sign",
                nac_convention=convention,
                migration_attestation=wrong_sign,
            ),
            (ValueError,),
            "quarantined",
        ),
        "corrected_pyscf_identity_is_full_overlap": bool(
            not full_pyscf.use_etfs
            and not full_pyscf.mult_ediff
            and "state-i-j" in full_pyscf.source_mapping_id
        ),
        "etf_identity_is_distinct": (
            full_pyscf.fingerprint() != etf_pyscf.fingerprint()
        ),
        "missing_snapshot_nac_identity_rejected": _raises_v233(
            lambda: require_snapshot_nac_identity_v233({}, full_pyscf),
            (ValueError,),
            "legacy data are quarantined",
        ),
        "mismatched_snapshot_nac_identity_rejected": _raises_v233(
            lambda: require_snapshot_nac_identity_v233(
                {"v233_nac_convention_fingerprint": etf_pyscf.fingerprint()},
                full_pyscf,
            ),
            (ValueError,),
            "identity mismatch",
        ),
        "direct_v233_capture_avoids_legacy_claim": bool(
            direct.migration_attestation is None
            and direct.nac_convention.fingerprint() == convention.fingerprint()
        ),
        "provider_numerical_identity_binds_all_conventions": bool(
            require_provider_numerical_identity_v233(provider, identity) == identity
            and identity.fingerprint() != changed_nac_identity.fingerprint()
            and identity.fingerprint() != changed_policy.fingerprint()
            and _raises_v233(
                lambda: require_provider_numerical_identity_v233(
                    LegacyProvider(), identity
                ),
                (ValueError,),
                "legacy providers are quarantined",
            )
        ),
        "dataset_fingerprint": first.dataset_fingerprint,
        "legacy_dataset_fingerprint": legacy.dataset_fingerprint,
        "nac_convention_fingerprint": convention.fingerprint(),
        "provider_numerical_identity_fingerprint": identity.fingerprint(),
    }


def _manifold_and_soc_controls_v233():
    even = AnalyticSingletTripletSOCProviderV220()
    odd = AnalyticDoubletSOCProviderV220()
    identity = np.eye(4, dtype=complex)
    even_report = audit_complete_manifold_transport_v233(
        identity, even.provenance.model_space, even.soc_symmetry_contract
    )
    odd_report = audit_complete_manifold_transport_v233(
        identity, odd.provenance.model_space, odd.soc_symmetry_contract
    )
    left_gauge = _block_unitary_v233(23331)
    right_gauge = _block_unitary_v233(23333)
    gauged_report = audit_complete_manifold_transport_v233(
        left_gauge.conj().T @ right_gauge,
        odd.provenance.model_space,
        _transform_symmetry_v233(odd.soc_symmetry_contract, left_gauge),
        right_symmetry_contract=_transform_symmetry_v233(
            odd.soc_symmetry_contract, right_gauge
        ),
    )
    broken_kramers = audit_complete_manifold_transport_v233(
        np.diag([1.0, 0.8, 1.0, 1.0]),
        odd.provenance.model_space,
        odd.soc_symmetry_contract,
    )
    swapped = identity.copy()
    swapped[:, [0, 1]] = swapped[:, [1, 0]]
    leakage = audit_complete_manifold_transport_v233(
        swapped, even.provenance.model_space, even.soc_symmetry_contract
    )
    incomplete_contract = SOCSymmetryContractV221(
        electron_parity="odd",
        time_reversal_matrix=odd.time_reversal_matrix,
        projectors={"doublet_1": odd.projectors["doublet_1"]},
    )
    incomplete = audit_complete_manifold_transport_v233(
        identity, odd.provenance.model_space, incomplete_contract
    )
    even_convention = analytic_soc_convention_v233(even)
    odd_convention = analytic_soc_convention_v233(odd)
    even_soc = audit_molecular_soc_convention_v233(
        even.components(np.asarray([0.17])),
        even.provenance,
        even.soc_symmetry_contract,
        even_convention,
    )
    odd_soc = audit_molecular_soc_convention_v233(
        odd.components(np.asarray([-0.11])),
        odd.provenance,
        odd.soc_symmetry_contract,
        odd_convention,
    )
    wrong_order = replace(
        even_convention, state_order=tuple(reversed(even_convention.state_order))
    )
    wrong_order_report = audit_molecular_soc_convention_v233(
        even.components(np.asarray([0.17])),
        even.provenance,
        even.soc_symmetry_contract,
        wrong_order,
    )
    wrong_prefactor = replace(
        even_convention, prefactor_convention="untrusted alternative prefactor"
    )
    return {
        "even_complete_manifolds": even_report.passed,
        "odd_complete_kramers_manifolds": odd_report.passed,
        "arbitrary_endpoint_doublet_gauges": gauged_report.passed,
        "broken_kramers_overlap_detected": bool(
            not broken_kramers.checks["time_reversal_covariance"]
        ),
        "competing_manifold_leakage_detected": bool(
            not leakage.checks["assigned_manifold_retention"]
            and not leakage.checks["competing_manifold_leakage"]
        ),
        "incomplete_projector_family_detected": bool(
            not incomplete.checks["left_complete_symmetry_contract"]
            and not incomplete.checks["right_complete_symmetry_contract"]
        ),
        "even_soc_matrix_convention": even_soc.passed,
        "odd_soc_matrix_convention": odd_soc.passed,
        "soc_state_order_mismatch_detected": bool(
            not wrong_order_report.checks["exact_state_order"]
        ),
        "soc_prefactor_trust_anchor_enforced": _raises_v233(
            lambda: require_exact_molecular_soc_convention_v233(
                wrong_prefactor, even_convention
            ),
            (ValueError,),
            "differs from the trusted",
        ),
        "even_report": even_report.as_dict(),
        "odd_report": odd_report.as_dict(),
        "gauged_doublet_report": gauged_report.as_dict(),
        "even_soc_convention_fingerprint": even_convention.fingerprint(),
        "odd_soc_convention_fingerprint": odd_convention.fingerprint(),
    }


def _runtime_profile_controls_v233(runtime):
    locked_profile = release_locked_runtime_profile_v233()
    compatible_profile = scientifically_compatible_runtime_profile_v233()
    locked = assess_runtime_compatibility_v233(runtime, locked_profile)
    compatible = assess_runtime_compatibility_v233(runtime, compatible_profile)
    portable = dict(runtime)
    portable["python_executable_sha256"] = "0" * 64
    portable["platform"] = "Linux-different-kernel"
    changed_locked = assess_runtime_compatibility_v233(portable, locked_profile)
    changed_compatible = assess_runtime_compatibility_v233(
        portable, compatible_profile
    )
    return {
        "release_locked_profile_passes_canonical_runtime": locked.compatible,
        "scientific_profile_passes_canonical_runtime": compatible.compatible,
        "portable_runtime_does_not_claim_byte_identity": bool(
            not changed_locked.compatible and changed_compatible.compatible
        ),
        "profile_fingerprints_are_distinct": (
            locked.profile_fingerprint != compatible.profile_fingerprint
        ),
        "release_locked_report": locked.as_dict(),
        "scientifically_compatible_report": compatible.as_dict(),
    }


def run_v0233_release_benchmark(thresholds=V233AcceptanceThresholds()):
    inherited = run_v0232_release_benchmark()
    with tempfile.TemporaryDirectory(prefix="gnd-v233-release-") as directory:
        transport = _transport_controls_v233()
        replay_nac = _replay_and_nac_controls_v233(Path(directory) / "replay")
    manifold_soc = _manifold_and_soc_controls_v233()
    runtime_profiles = _runtime_profile_controls_v233(
        inherited["pyscf_runtime_evidence"]["runtime"]
    )

    new_checks = {
        **{
            f"transport::{name}": bool(value)
            for name, value in transport.items()
            if name not in {"diagnostics", "low_retention_diagnostics"}
        },
        **{
            f"replay_nac::{name}": bool(value)
            for name, value in replay_nac.items()
            if name
            not in {
                "dataset_fingerprint",
                "legacy_dataset_fingerprint",
                "nac_convention_fingerprint",
                "provider_numerical_identity_fingerprint",
                # Reported as supporting evidence; the two byte-level
                # determinism gates already bind this derived value.
                "migration_fingerprint_deterministic",
            }
        },
        **{
            f"manifold_soc::{name}": bool(value)
            for name, value in manifold_soc.items()
            if name
            not in {
                "even_report",
                "odd_report",
                "gauged_doublet_report",
                "even_soc_convention_fingerprint",
                "odd_soc_convention_fingerprint",
            }
        },
        **{
            f"runtime_profile::{name}": bool(value)
            for name, value in runtime_profiles.items()
            if name
            not in {
                "release_locked_report",
                "scientifically_compatible_report",
            }
        },
    }
    if len(new_checks) != thresholds.expected_new_gates:
        raise AssertionError(
            f"v0.23.3 campaign defines {len(new_checks)} rather than exactly "
            f"{thresholds.expected_new_gates} new gates."
        )
    inherited_checks = {
        f"inherited_v0232::{name}": bool(value)
        for name, value in inherited["acceptance"]["checks"].items()
    }
    if len(inherited_checks) != thresholds.expected_inherited_gates:
        raise AssertionError("v0.23.3 must inherit exactly 168 v0.23.2 gates.")
    checks = {**inherited_checks, **new_checks}
    if len(checks) != thresholds.expected_total_gates:
        raise AssertionError("v0.23.3 campaign must define exactly 208 total gates.")
    return {
        "release": "v0.23.3",
        "theme": (
            "finite-manifold overlap/transport separation, replay format 2, "
            "legacy NAC quarantine, complete multiplet transport, frozen SOC "
            "matrix conventions, and separate runtime profiles"
        ),
        "transport_controls": transport,
        "replay_and_nac_controls": replay_nac,
        "manifold_and_soc_controls": manifold_soc,
        "runtime_profile_controls": runtime_profiles,
        "claims": {
            "finite_manifold_unitary_transport_validated": True,
            "trajectory_overlap_quality_policy_validated": True,
            "replay_format_two_and_migration_validated": True,
            "legacy_NAC_data_quarantine_validated": True,
            "complete_even_and_odd_manifold_transport_validated": True,
            "molecular_SOC_matrix_convention_frozen": True,
            "runtime_identity_and_compatibility_profiles_separated": True,
            "real_PySCF_spin_free_runtime_inherited": True,
            "physical_analytic_SOC_inherited": True,
            "external_molecular_SOC_snapshot_admitted": False,
            "live_molecular_SOC_backend_admitted": False,
            "ab_initio_SOC_validated": False,
            "live_PySCF_SOC_runtime_validated": False,
        },
        "inherited_v0232": inherited,
        "acceptance": {
            "passed": bool(all(checks.values())),
            "checks": checks,
            "inherited_gate_count": len(inherited_checks),
            "new_gate_count": len(new_checks),
            "total_gate_count": len(checks),
            "thresholds": asdict(thresholds),
        },
    }
