import numpy as np

from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from gaussian_dynamics.curved_basis_adaptation_v292 import (
    ConstantConnectionV292, CurvedBasisSettingsV292, CurvedTDVPSnapshotV292,
    SVDReceiptV292, actual_linear_residual_reduction_v292,
    candidate_coefficient_block_v292, evaluate_curved_spawn_candidate_v292,
    finite_difference_full_jacobian_v292, no_event_reduction_v292,
    spawn_curved_candidate_v292, tensor_legendre_quadrature_v292,
)


def fixture():
    state = CorrelatedGaussianSpinorStateV270(
        q=np.array([[-0.35, 0.20]]), p=np.array([[0.18, -0.11]]),
        width_matrices=np.array([[[1.10, 0.18],[0.18,0.82]]]),
        chirp_matrices=np.array([[[0.07,-0.03],[-0.03,0.04]]]),
        coefficients=np.array([[0.82+0.10j, 0.42-0.16j]]),
    ).normalized()
    sy=np.array([[0,-1j],[1j,0]],complex); sz=np.array([[1,0],[0,-1]],complex)
    conn=ConstantConnectionV292(np.array([0.22j*sy, 0.17j*sz])).validate()
    quad=tensor_legendre_quadrature_v292(((-4,4),(-4,4)), (18,18))
    return state,conn,quad


def synthetic_snapshot(state,conn,quad,candidate):
    j=finite_difference_full_jacobian_v292(state,conn,quad,relative_step=3e-6)
    k=candidate_coefficient_block_v292(candidate,state.nstate,conn,quad)
    # Deliberately create a defect with a component in a genuinely new candidate direction.
    residual=0.16*k[:,0]-0.11*k[:,1]+0.03j*k[:,0]
    # A true McLachlan defect is real-orthogonal to every current tangent:
    # Re(J^dagger r)=0.  Enforce that defining condition in the synthetic oracle.
    jr=np.vstack([j.real,j.imag])
    rr=np.concatenate([residual.real,residual.imag])
    correction=np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]
    rr=rr+jr@correction
    n=j.shape[0]
    residual=rr[:n]+1j*rr[n:]
    m=np.real(j.conj().T@j)
    rhs=np.zeros(j.shape[1])
    s=np.linalg.svd(m,compute_uv=False)
    rec=SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0.0,0.0)
    return CurvedTDVPSnapshotV292(j,m,rhs,np.zeros(j.shape[1]),residual,rec).validate()


def test_predicted_gain_matches_independent_augmented_least_squares():
    state,conn,quad=fixture()
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    snap=synthetic_snapshot(state,conn,quad,candidate)
    settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
    ev=evaluate_curved_spawn_candidate_v292(state,conn,quad,snap,candidate,settings=settings)
    k=candidate_coefficient_block_v292(candidate,state.nstate,conn,quad)
    actual=actual_linear_residual_reduction_v292(snap,k)
    assert ev.admitted
    assert ev.predicted_residual_reduction > 0
    assert abs(ev.predicted_residual_reduction-actual)/actual < 2e-6


def test_existing_packet_coefficient_block_adds_no_new_tangent_space():
    state,conn,quad=fixture()
    # candidate exactly duplicates current geometry
    from gaussian_dynamics.correlated_basis_adaptation_v270 import MultidimensionalSpawnCandidateV270
    c=MultidimensionalSpawnCandidateV270(state.q[0],state.p[0],state.width_matrices[0],state.chirp_matrices[0],0,"position",0,1)
    snap=synthetic_snapshot(state,conn,quad,generate_multidimensional_spawn_candidates_v270(state)[0])
    ev=evaluate_curved_spawn_candidate_v292(state,conn,quad,snap,c,settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12))
    assert not ev.admitted
    assert ev.predicted_residual_reduction < 1e-8


def test_spawn_projection_preserves_curved_wavefunction_to_quadrature_precision():
    state,conn,quad=fixture()
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    out,receipt=spawn_curved_candidate_v292(state,candidate,conn,quad)
    assert out.ngaussian==2
    assert receipt.relative_loss < 5e-20
    assert receipt.fidelity > 1-5e-13
    assert np.linalg.norm(out.coefficients[-1]) < 1e-10


def test_no_event_is_exact_object_identity():
    state,_,_=fixture()
    assert no_event_reduction_v292(state) is state


def test_constant_connection_is_antihermitian_and_wilson_is_unitary():
    _,conn,_=fixture()
    l=conn.wilson(np.array([0.4,-0.3]))
    assert np.linalg.norm(l.conj().T@l-np.eye(2)) < 2e-13


def test_degenerate_width_candidate_generation_fails_closed():
    state,_,_=fixture()
    iso=CorrelatedGaussianSpinorStateV270(state.q,state.p,np.array([np.eye(2)]),state.chirp_matrices,state.coefficients).normalized()
    assert generate_multidimensional_spawn_candidates_v270(iso)==()


def identity_hamiltonian(state, psi, quadrature):
    return np.asarray(psi, dtype=complex)


def dummy_snapshot_for_state(state,conn,quad):
    j=finite_difference_full_jacobian_v292(state,conn,quad,relative_step=3e-6)
    m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False)
    return CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),np.zeros(j.shape[0],complex),SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0,0)).validate()


def test_merge_uses_full_curved_projection_and_energy_guard():
    from gaussian_dynamics.curved_basis_adaptation_v292 import adapt_curved_basis_once_v292
    state,conn,quad=fixture()
    # Two identical packet subspaces with split amplitudes are exactly mergeable.
    two=CorrelatedGaussianSpinorStateV270(
        np.vstack([state.q,state.q]),np.vstack([state.p,state.p]),
        np.concatenate([state.width_matrices,state.width_matrices]),np.concatenate([state.chirp_matrices,state.chirp_matrices]),
        np.vstack([0.7*state.coefficients,0.3*state.coefficients]),0.0).normalized()
    snap=dummy_snapshot_for_state(two,conn,quad)
    ev=adapt_curved_basis_once_v292(two,conn,quad,snap,identity_hamiltonian,packet_ages=(100,100))
    assert ev.event_kind=="merge"
    assert ev.after.ngaussian==1
    assert ev.projection.relative_loss < 1e-18
    assert abs(ev.energy_jump_hartree) < 1e-12


def test_prune_requires_age_population_projection_and_energy_gates():
    from gaussian_dynamics.curved_basis_adaptation_v292 import adapt_curved_basis_once_v292
    state,conn,quad=fixture()
    q2=np.array([[2.8,-2.6]])
    two=CorrelatedGaussianSpinorStateV270(
        np.vstack([state.q,q2]),np.vstack([state.p,[[0.0,0.0]]]),
        np.concatenate([state.width_matrices,state.width_matrices]),np.concatenate([state.chirp_matrices,state.chirp_matrices]),
        np.vstack([state.coefficients,[[1e-7+0j,0j]]]),0.0).normalized()
    snap=dummy_snapshot_for_state(two,conn,quad)
    settings=CurvedBasisSettingsV292(inherited=__import__('gaussian_dynamics.correlated_basis_adaptation_v270',fromlist=['ControlledMultidimensionalBasisSettingsV270']).ControlledMultidimensionalBasisSettingsV270(prune_population_threshold=1e-10,minimum_prune_age=2,minimum_merge_overlap=.999999))
    ev=adapt_curved_basis_once_v292(two,conn,quad,snap,identity_hamiltonian,packet_ages=(10,10),settings=settings)
    assert ev.event_kind=="prune"
    assert ev.after.ngaussian==1
    assert ev.projection.relative_loss < settings.inherited.maximum_projection_loss


def test_spawn_event_marks_new_shape_dormant():
    from gaussian_dynamics.curved_basis_adaptation_v292 import adapt_curved_basis_once_v292
    state,conn,quad=fixture()
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    snap=synthetic_snapshot(state,conn,quad,candidate)
    settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
    ev=adapt_curved_basis_once_v292(state,conn,quad,snap,identity_hamiltonian,settings=settings)
    assert ev.event_kind=="spawn"
    assert ev.after.ngaussian==2
    assert ev.active_shape_mask_after.tolist()==[True,False]
    assert np.linalg.norm(ev.after.coefficients[-1]) < 1e-10
