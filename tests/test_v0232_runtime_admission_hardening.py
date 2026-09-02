from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gaussian_dynamics.molecular_soc_admission_v232 import (
    audit_molecular_soc_provider_v232,
)
from gaussian_dynamics.molecular_soc_dossier_v231 import (
    BackendRuntimeAttestationV231,
    write_molecular_soc_dossier_v231,
    write_raw_json_artifact_v231,
)
from gaussian_dynamics.molecular_soc_runtime_v232 import (
    BackendAdmissionPolicyV232,
    BackendArtifactValidationProofV232,
    BackendMethodIdentityV232,
    CONVERGENCE_METADATA_KEY_V232,
    ConvergenceMetadataV232,
    ReceiptExecutionEvidenceV232,
    RuntimeProbeRecordV232,
    convergence_from_snapshot_v232,
)
from gaussian_dynamics.pyscf_soc_adapter_v232 import (
    PYSCF_NAC_CONVENTION_V232,
    PySCFMethodSpecificCapabilitiesV232,
    validate_pyscf_engine_contract_v232,
)
from gaussian_dynamics.v231_benchmark import build_v231_admission_bundle


def _pyscf_identity(**changes):
    fields = dict(
        backend_name="PySCF",
        backend_version="2.13.1",
        source_kind="live_ab_initio",
        adapter_name="trusted-pyscf-soc",
        adapter_version="1",
        electronic_method="SA-CASSCF/SI-SOC",
        basis="cc-pVDZ",
        active_space="CAS(2,2)",
        soc_operator="Breit-Pauli SOMF",
        scalar_relativistic_method="none",
        derivative_method="analytic spin-free + physical SOC",
        nac_convention=PYSCF_NAC_CONVENTION_V232,
    )
    fields.update(changes)
    return BackendMethodIdentityV232(**fields).validate()


def test_engine_structure_is_checked_before_capability_declarations():
    class DeclarationTrap:
        method_identity = _pyscf_identity()

        @property
        def capabilities(self):
            raise AssertionError("capabilities were consulted before structure")

    with pytest.raises(TypeError, match="incomplete"):
        validate_pyscf_engine_contract_v232(
            DeclarationTrap(),
            runtime_version="2.13.1",
            expected_identity=_pyscf_identity(),
        )


def test_engine_rejects_noncallable_methods_and_exact_method_mismatch():
    class Engine:
        method_identity = _pyscf_identity()
        capabilities = PySCFMethodSpecificCapabilitiesV232(
            True, True, True, True, True, True, True
        )
        components = None

        def evaluate_snapshot(self):
            pass

        def snapshot_overlap(self):
            pass

        def write_raw_artifacts(self):
            pass

        def validate_raw_artifacts_v232(self):
            pass

    with pytest.raises(TypeError, match="noncallable components"):
        validate_pyscf_engine_contract_v232(
            Engine(), runtime_version="2.13.1", expected_identity=_pyscf_identity()
        )

    Engine.components = lambda self: None
    with pytest.raises(ValueError, match="trusted identity"):
        validate_pyscf_engine_contract_v232(
            Engine(),
            runtime_version="2.13.1",
            expected_identity=_pyscf_identity(basis="aug-cc-pVDZ"),
        )


def test_convergence_vocabulary_is_explicit_and_has_no_override_merge():
    legacy = SimpleNamespace(
        point=SimpleNamespace(metadata={"scf_converged": True}),
        metadata={"derivatives_converged": True},
    )
    with pytest.raises(ValueError, match="lacks canonical"):
        convergence_from_snapshot_v232(legacy)

    convergence = ConvergenceMetadataV232(
        **{name: np.bool_(True) for name in (
            "scf",
            "correlated_wavefunction",
            "state_interaction_soc",
            "spin_free_gradients",
            "soc_derivatives",
            "derivative_connections",
            "many_electron_overlaps",
        )}
    )
    snapshot = SimpleNamespace(
        point=SimpleNamespace(metadata={}),
        metadata={CONVERGENCE_METADATA_KEY_V232: convergence.as_dict()},
    )
    assert convergence_from_snapshot_v232(snapshot).complete

    snapshot.point.metadata[CONVERGENCE_METADATA_KEY_V232] = convergence.as_dict()
    with pytest.raises(ValueError, match="duplicate"):
        convergence_from_snapshot_v232(snapshot)


def _runtime_identity_from_bundle(bundle):
    identity = bundle["provider"].molecular_soc_contract.identity
    return BackendMethodIdentityV232(
        backend_name=identity.backend_name,
        backend_version=identity.backend_version,
        source_kind=identity.source_kind,
        adapter_name="fixture-parser-adapter",
        adapter_version="1",
        electronic_method=identity.electronic_method,
        basis=identity.basis,
        active_space=identity.active_space,
        soc_operator=identity.soc_operator,
        scalar_relativistic_method=identity.scalar_relativistic_method,
        derivative_method=identity.derivative_method,
        nac_convention=PYSCF_NAC_CONVENTION_V232,
    ).validate()


def _attested_bundle(
    tmp_path, *, source_kind="live_ab_initio", fresh_execution=True
):
    bundle = build_v231_admission_bundle(
        tmp_path / source_kind, source_kind=source_kind
    )
    identity = bundle["provider"].molecular_soc_contract.identity
    method_identity = _runtime_identity_from_bundle(bundle)
    runtime_probe = RuntimeProbeRecordV232(
        method_identity=method_identity,
        environment_sha256=identity.environment_sha256,
        calculation_input_sha256=identity.calculation_input_sha256,
        replay_dataset_fingerprint=bundle["dataset"].dataset_fingerprint,
        runtime_imported=True,
    ).validate()
    probe_artifact = write_raw_json_artifact_v231(
        bundle["directory"],
        name="runtime_probe_v232",
        relative_path="raw/runtime_probe_v232.json",
        role="runtime_probe",
        payload=runtime_probe.as_dict(),
    )
    attestation = BackendRuntimeAttestationV231(
        runtime_name=method_identity.backend_name,
        runtime_version=method_identity.backend_version,
        adapter_name=method_identity.adapter_name,
        adapter_version=method_identity.adapter_version,
        environment_sha256=identity.environment_sha256,
        runtime_probe_artifact=probe_artifact.name,
        runtime_imported=True,
        method_specific_soc_implemented=True,
        soc_derivatives_implemented=True,
        wavefunction_overlaps_implemented=True,
        artifact_parser_validated=True,
        fresh_execution_observed=fresh_execution,
    ).validate()
    dossier = replace(
        bundle["dossier"],
        artifacts=(*bundle["dossier"].artifacts, probe_artifact),
        runtime_attestation=attestation,
    )
    dossier = write_molecular_soc_dossier_v231(
        bundle["directory"],
        dossier,
        dataset=bundle["dataset"],
        identity=identity,
        overwrite=True,
    )
    bundle.update(
        dossier=dossier,
        method_identity=method_identity,
        probe_artifact=probe_artifact,
    )
    return bundle


class _TrustedFixtureValidator:
    parser_name = "strict-fixture-parser"
    parser_version = "1"

    def __init__(self, method_identity, *, alter=None, fresh_execution=True):
        self.method_identity = method_identity
        self.alter = alter
        self.fresh_execution = fresh_execution
        self.calls = 0

    def validate_raw_artifacts_v232(
        self, *, dossier, bundle_directory, dataset, execution_challenge
    ):
        self.calls += 1
        artifacts = {item.name: item for item in dossier.artifacts}
        convergence = ConvergenceMetadataV232(
            scf=True,
            correlated_wavefunction=True,
            state_interaction_soc=True,
            spin_free_gradients=True,
            soc_derivatives=True,
            derivative_connections=True,
            many_electron_overlaps=True,
        )
        proof = BackendArtifactValidationProofV232(
            method_identity=self.method_identity,
            parser_name=self.parser_name,
            parser_version=self.parser_version,
            execution_challenge=execution_challenge,
            replay_dataset_fingerprint=dataset.dataset_fingerprint,
            dossier_fingerprint=dossier.fingerprint(),
            environment_sha256=dossier.runtime_attestation.environment_sha256,
            calculation_input_sha256=artifacts[
                dossier.calculation_template_artifact
            ].sha256,
            runtime_probe_artifact=dossier.runtime_attestation.runtime_probe_artifact,
            runtime_probe_sha256=artifacts[
                dossier.runtime_attestation.runtime_probe_artifact
            ].sha256,
            parsed_output_artifacts=tuple(
                (receipt.output_artifact, artifacts[receipt.output_artifact].sha256)
                for receipt in dossier.receipts
            ),
            receipt_evidence=tuple(
                ReceiptExecutionEvidenceV232(
                    receipt.record_id, receipt.output_artifact, convergence
                )
                for receipt in dossier.receipts
            ),
            parser_executed=True,
            fresh_execution_observed=self.fresh_execution,
        )
        if self.alter == "bare_boolean":
            return True
        if self.alter == "wrong_challenge":
            proof = replace(proof, execution_challenge="f" * 64)
        if self.alter == "missing_output":
            proof = replace(
                proof, parsed_output_artifacts=proof.parsed_output_artifacts[:-1]
            )
        return proof.validate()


class _SelfDeclaredImpostor(_TrustedFixtureValidator):
    pass


def _policy(bundle):
    return BackendAdmissionPolicyV232(
        expected_identity=bundle["method_identity"],
        trusted_validator_type=_TrustedFixtureValidator,
        parser_name=_TrustedFixtureValidator.parser_name,
        parser_version=_TrustedFixtureValidator.parser_version,
    ).validate()


def _audit(bundle, validator, *, requirement="live"):
    return audit_molecular_soc_provider_v232(
        bundle["provider"],
        bundle["center"],
        bundle["dossier"],
        requirement=requirement,
        policy=_policy(bundle),
        backend_validator=validator,
        execution_challenge="a" * 64,
        bundle_directory=bundle["directory"],
    )


def test_typed_proof_is_bound_but_explicit_negative_control_remains_blocked(tmp_path):
    bundle = _attested_bundle(tmp_path)
    validator = _TrustedFixtureValidator(bundle["method_identity"])
    report = _audit(bundle, validator)

    assert not report.live_backend_admitted
    assert not report.passed
    assert validator.calls == 1
    assert report.checks["component_resolved_convergence"]
    assert report.checks["exact_parsed_output_inventory"]
    assert not report.checks["explicit_negative_control_absent"]
    assert not report.checks["explicit_false_molecular_claim_absent"]


@pytest.mark.parametrize("alter", ["bare_boolean", "wrong_challenge", "missing_output"])
def test_boolean_and_incompletely_bound_parser_results_fail_closed(tmp_path, alter):
    bundle = _attested_bundle(tmp_path)
    report = _audit(
        bundle, _TrustedFixtureValidator(bundle["method_identity"], alter=alter)
    )
    assert not report.live_backend_admitted
    assert not report.passed


def test_self_declared_validator_cannot_select_its_own_trust_type(tmp_path):
    bundle = _attested_bundle(tmp_path)
    impostor = _SelfDeclaredImpostor(bundle["method_identity"])
    report = _audit(bundle, impostor)

    assert not report.live_backend_admitted
    assert not report.checks["trusted_validator_exact_type"]
    assert impostor.calls == 0


def test_live_attestation_cannot_disagree_with_fresh_execution_proof(tmp_path):
    bundle = _attested_bundle(tmp_path)
    attestation = replace(
        bundle["dossier"].runtime_attestation, fresh_execution_observed=False
    )
    bundle["dossier"] = replace(
        bundle["dossier"], runtime_attestation=attestation
    )
    validator = _TrustedFixtureValidator(bundle["method_identity"])
    report = _audit(bundle, validator)

    assert not report.live_backend_admitted
    assert not report.checks["runtime_attestation_live_execution"]
    assert validator.calls == 0


def test_external_parser_executes_without_fresh_live_execution(tmp_path):
    bundle = _attested_bundle(
        tmp_path,
        source_kind="external_ab_initio_snapshot",
        fresh_execution=False,
    )
    validator = _TrustedFixtureValidator(
        bundle["method_identity"], fresh_execution=False
    )
    report = _audit(bundle, validator, requirement="external")

    assert validator.calls == 1
    assert report.checks["typed_parser_execution_proof"]
    assert report.checks["parser_actually_invoked"]
    assert not report.checks["runtime_attestation_live_execution"]
    assert not report.checks["fresh_backend_execution_observed"]
    # This remains a labelled synthetic negative control, never a real admission.
    assert not report.external_snapshot_admitted


def test_live_parser_does_not_run_without_fresh_execution_attestation(tmp_path):
    bundle = _attested_bundle(tmp_path, fresh_execution=False)
    validator = _TrustedFixtureValidator(
        bundle["method_identity"], fresh_execution=False
    )
    report = _audit(bundle, validator, requirement="live")

    assert validator.calls == 0
    assert not report.checks["runtime_attestation_live_execution"]
    assert not report.live_backend_admitted
