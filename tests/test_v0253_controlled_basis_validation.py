from dataclasses import replace

import pytest

from gaussian_dynamics.controlled_basis_adaptation_v253 import (
    V253_CONTROLLED_BASIS_CLAIMS,
)
from gaussian_dynamics.controlled_basis_validation_v253 import (
    CONTROLLED_BASIS_VALIDATION_SCHEMA_V253,
    run_controlled_basis_validation_evidence_v253,
)


@pytest.fixture(scope="module")
def evidence():
    return run_controlled_basis_validation_evidence_v253()


def test_validation_has_exactly_sixty_boolean_passing_gates(evidence):
    assert len(evidence.audit.checks) == 60
    assert all(type(value) is bool for value in evidence.audit.checks.values())
    assert all(evidence.audit.checks.values())
    assert evidence.audit.passed is True


def test_validation_schema_decisions_and_claim_boundary_are_explicit(evidence):
    payload = evidence.as_dict()
    assert payload["schema"] == CONTROLLED_BASIS_VALIDATION_SCHEMA_V253
    assert "dPsi/dt+iHPsi" in payload["decisions"]["residual_score"]
    assert "full-SVD" in payload["decisions"]["projection_policy"]
    assert V253_CONTROLLED_BASIS_CLAIMS["controlled_residual_driven_spawning_validated"] is True
    assert V253_CONTROLLED_BASIS_CLAIMS["general_aims_branching_validated"] is False
    assert V253_CONTROLLED_BASIS_CLAIMS["real_pyscf_soc_trajectory_admitted"] is False


def test_evidence_binds_spawn_activation_prune_merge_and_reduction(evidence):
    assert evidence.odd_spawn_event.event_kind == "spawn"
    assert evidence.activation_step.active_shape_mask.tolist() == [True, True, False]
    assert evidence.prune_event.event_kind == "prune"
    assert evidence.merge_event.event_kind == "merge"
    assert evidence.controlled_trajectory.event_counts["spawn"] == 1
    assert evidence.no_event_trajectory.event_counts["none"] == 1


def test_evidence_fingerprint_is_deterministic(evidence):
    first = evidence.fingerprint()
    second = run_controlled_basis_validation_evidence_v253().fingerprint()
    assert first == second
    assert len(first) == 64
    assert all(character in "0123456789abcdef" for character in first)


def test_tampered_validation_audit_is_rejected(evidence):
    checks = dict(evidence.audit.checks)
    checks["odd_spawn_event_accepted"] = False
    with pytest.raises(ValueError, match="audit result is inconsistent"):
        replace(evidence.audit, checks=checks).validate()
    with pytest.raises(ValueError, match="requires exactly 60"):
        replace(evidence.audit, checks={}).validate()
