import pytest
from gaussian_dynamics.moving_frame_evidence_v280 import MOVING_FRAME_EVIDENCE_SCHEMA_V280,run_moving_frame_evidence_v280
@pytest.fixture(scope='module')
def evidence(): return run_moving_frame_evidence_v280()
def test_v0280_evidence_passes_all_frozen_gates(evidence): assert evidence.passed is True; assert evidence.check_count==50; assert all(evidence.checks.values())
def test_v0280_evidence_schema_and_fingerprints(evidence): assert MOVING_FRAME_EVIDENCE_SCHEMA_V280.endswith('v0.28.0'); assert len(evidence.trajectory_fingerprint)==64; assert len(evidence.lattice_fingerprint)==64; assert len(evidence.fingerprint())==64
def test_v0280_evidence_independent_lattice_defects_are_machine_precision(evidence):
    m=evidence.metrics; assert m['lattice_similarity_residual']<2e-15; assert m['lattice_action_covariance_residual']<2e-15; assert m['lattice_propagation_covariance_residual']<3e-15
def test_v0280_evidence_claim_boundaries_are_explicit(evidence):
    c=evidence.claims; assert c['validated']['flat_coordinate_dependent_gauge_covariance'] is True; assert c['not_validated']['nonzero_curvature_connections'] is True; assert c['not_validated']['live_molecular_soc_trajectories'] is True; assert c['not_validated']['full_aims_branching_semantics'] is True
