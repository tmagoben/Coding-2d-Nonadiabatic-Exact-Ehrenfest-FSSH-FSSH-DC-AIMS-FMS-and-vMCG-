from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics import (
    AnalyticSingletTripletSOCProviderV220,
    MolecularSOCAdmissionContractV230,
    MolecularSOCBackendIdentityV230,
    MolecularSOCCapabilitiesV230,
    MolecularSOCValidationEvidenceV230,
    provenance_with_molecular_soc_contract_v230,
    require_trajectory_ready_molecular_soc_v230,
)
from gaussian_dynamics.v230_benchmark import _fixture_contract_v230


def test_capability_tiers_are_derived_from_complete_requirements():
    static = MolecularSOCCapabilitiesV230(static_soc=True)
    trajectory = MolecularSOCCapabilitiesV230(
        static_soc=True,
        spin_free_derivatives=True,
        soc_derivatives=True,
        derivative_connections=True,
        cross_geometry_overlaps=True,
    )

    assert static.tier == "static_soc"
    assert not static.trajectory_ready
    assert trajectory.tier == "trajectory_ready"
    assert trajectory.trajectory_ready


def test_analytic_derivative_declaration_requires_soc_derivatives():
    with pytest.raises(ValueError, match="requires the SOC-derivative"):
        MolecularSOCCapabilitiesV230(
            static_soc=True, analytic_soc_derivatives=True
        ).validate()


def test_backend_identity_freezes_atomic_units_and_electron_parity():
    identity = _fixture_contract_v230().identity

    assert identity.electron_parity == "even"
    with pytest.raises(ValueError, match="bohr"):
        replace(identity, geometry_unit="angstrom").validate()
    with pytest.raises(ValueError, match="positive integer"):
        replace(identity, electron_count=0).validate()


def test_evidence_cannot_be_partially_declared():
    with pytest.raises(ValueError, match="incomplete"):
        MolecularSOCValidationEvidenceV230(
            independent_reference_id="reference-only"
        ).validate()
    with pytest.raises(ValueError, match="one value"):
        MolecularSOCValidationEvidenceV230(
            basis_levels=("a", "b", "c"),
            basis_changes=(1.0e-4,),
            basis_tolerance=1.0e-5,
        ).validate()


def test_real_admission_readiness_requires_source_and_all_evidence():
    reference = MolecularSOCValidationEvidenceV230(
        independent_reference_id="independent data DOI",
        independent_reference_error=1.0e-6,
        independent_reference_tolerance=1.0e-5,
        basis_levels=("small", "large"),
        basis_changes=(1.0e-6,),
        basis_tolerance=1.0e-5,
        method_levels=("method-a", "method-b"),
        method_changes=(2.0e-6,),
        method_tolerance=1.0e-5,
        translation_residual=1.0e-9,
        rotation_residual=2.0e-9,
        frame_invariance_tolerance=1.0e-8,
        tracking_minimum_overlap=0.92,
        tracking_minimum_margin=0.18,
        tracking_overlap_threshold=0.80,
        tracking_margin_threshold=0.10,
    )
    fixture = _fixture_contract_v230(evidence=reference)
    external = _fixture_contract_v230(
        source_kind="external_ab_initio_snapshot",
        evidence=reference,
    )

    assert not fixture.real_backend_admission_ready
    assert external.real_backend_admission_ready


def test_real_source_requires_explicit_nuclear_and_environment_identity():
    fixture = _fixture_contract_v230()

    with pytest.raises(ValueError, match="nuclear identity"):
        replace(
            fixture.identity,
            source_kind="external_ab_initio_snapshot",
        ).validate()


def test_electron_count_must_match_soc_symmetry_parity():
    provider = AnalyticSingletTripletSOCProviderV220()
    contract = _fixture_contract_v230()
    odd_identity = replace(contract.identity, electron_count=3)

    with pytest.raises(ValueError, match="electron count"):
        replace(contract, identity=odd_identity).validate(
            provider.soc_symmetry_contract
        )


def test_molecular_contract_is_part_of_operator_provenance_identity():
    provider = AnalyticSingletTripletSOCProviderV220()
    contract = _fixture_contract_v230()
    changed = replace(
        contract,
        coordinate_definition="a different generalized coordinate definition",
    )
    first = provenance_with_molecular_soc_contract_v230(
        provider.provenance, contract
    )
    second = provenance_with_molecular_soc_contract_v230(
        provider.provenance, changed
    )

    assert first.fingerprint() != second.fingerprint()
    assert first.parameters["v230_molecular_soc_contract_fingerprint"] == (
        contract.fingerprint()
    )


def test_static_only_provider_is_rejected_for_moving_nuclei():
    provider = AnalyticSingletTripletSOCProviderV220()

    class StaticProbe:
        soc_symmetry_contract = provider.soc_symmetry_contract
        molecular_soc_contract = _fixture_contract_v230(
            capabilities=MolecularSOCCapabilitiesV230(static_soc=True)
        )

    with pytest.raises(ValueError, match="moving-nuclear dynamics"):
        require_trajectory_ready_molecular_soc_v230(StaticProbe())
