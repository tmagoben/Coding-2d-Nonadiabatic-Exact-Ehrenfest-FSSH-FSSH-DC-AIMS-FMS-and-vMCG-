import numpy as np
from scipy.linalg import expm

from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.correlated_basis_adaptation_v270 import generate_multidimensional_spawn_candidates_v270
from gaussian_dynamics.curved_basis_adaptation_v292 import (
    ConstantConnectionV292, CurvedBasisSettingsV292, CurvedTDVPSnapshotV292,
    SVDReceiptV292, actual_linear_residual_reduction_v292,
    coefficient_basis_v292, tensor_legendre_quadrature_v292,
)
from gaussian_dynamics.curved_geometry_v293 import (
    ConnectionFieldV293, CoordinateGaugeV293, PatchTransitionV293, PolylinePathV293,
    constant_connection_field_v293, gauge_transform_connection_v293,
    straight_path_v293, dogleg_path_v293, wilson_path_v293,
    plaquette_holonomy_v293, validate_patch_transition_v293,
    two_patch_transport_v293,
)
from gaussian_dynamics.curved_basis_adaptation_v293 import (
    coefficient_basis_v293, finite_difference_full_jacobian_v293,
    candidate_coefficient_block_v293, evaluate_curved_spawn_candidate_v293,
    spawn_curved_candidate_v293,
)


def paulis():
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.array([[1,0],[0,-1]],complex)
    return sx,sy,sz


def fixture_state():
    return CorrelatedGaussianSpinorStateV270(
        q=np.array([[-0.25,0.15]]), p=np.array([[0.12,-0.08]]),
        width_matrices=np.array([[[1.05,0.14],[0.14,0.86]]]),
        chirp_matrices=np.array([[[0.05,-0.02],[-0.02,0.03]]]),
        coefficients=np.array([[0.84+0.08j,0.43-0.14j]]),
    ).normalized()


def base_constant():
    _,sy,sz=paulis()
    gens=np.array([0.19j*sy,0.13j*sz])
    return gens,constant_connection_field_v293(gens)


def coordinate_gauge():
    _,_,sz=paulis(); x=0.35j*sz
    def theta(r): return 0.18*r[0]-0.11*r[1]+0.07*r[0]*r[1]
    def grad(r): return np.array([0.18+0.07*r[1],-0.11+0.07*r[0]])
    def mat(r): return expm(theta(r)*x)
    def deriv(r):
        g=mat(r); return np.asarray([grad(r)[mu]*x@g for mu in range(2)])
    return CoordinateGaugeV293(2,2,mat,deriv,"smooth-gauge")


def test_constant_field_reduces_to_v292_straight_sections():
    state=fixture_state(); gens,field=base_constant(); old=ConstantConnectionV292(gens).validate()
    quad=tensor_legendre_quadrature_v292(((-2.2,2.2),(-2.2,2.2)),(5,5))
    a=coefficient_basis_v292(state,old,quad); b=coefficient_basis_v293(state,field,quad)
    assert np.max(np.abs(a-b)) < 5e-13


def test_coordinate_dependent_gauge_covariance_of_wilson_line():
    _,source=base_constant(); gauge=coordinate_gauge(); target=gauge_transform_connection_v293(source,gauge)
    a=np.array([-0.4,0.2]); b=np.array([0.55,-0.35])
    w=wilson_path_v293(source,straight_path_v293(a,b))
    wp=wilson_path_v293(target,straight_path_v293(a,b))
    expected=gauge.at(b).conj().T@w@gauge.at(a)
    assert np.linalg.norm(wp-expected) < 3e-12
    assert np.linalg.norm(wp.conj().T@wp-np.eye(2)) < 3e-12


def variable_field():
    sx,sy,sz=paulis()
    # A_x and A_y are anti-Hermitian and both coordinate dependent.
    def eval_(r):
        x,y=r
        ax=0.14j*sy + 0.055j*y*sz
        ay=0.11j*sx - 0.045j*x*sz
        return np.array([ax,ay])
    return ConnectionFieldV293(2,2,eval_,"variable-nonabelian").validate([np.zeros(2),np.array([.3,-.2])])


def test_nonzero_curvature_makes_path_family_physically_distinct():
    field=variable_field(); a=np.array([-.35,-.25]); b=np.array([.45,.40])
    straight=wilson_path_v293(field,straight_path_v293(a,b))
    dogleg=wilson_path_v293(field,dogleg_path_v293(a,b,axis=0))
    assert np.linalg.norm(straight-dogleg) > 1e-4
    assert np.linalg.norm(straight.conj().T@straight-np.eye(2)) < 2e-10
    assert np.linalg.norm(dogleg.conj().T@dogleg-np.eye(2)) < 2e-10


def test_small_plaquette_recovers_local_curvature_with_sign():
    field=variable_field(); r=np.array([0.13,-0.21]); h=1.5e-3
    w=plaquette_holonomy_v293(field,r,0,1,h)
    f=field.curvature(r,0,1,step=2e-6)
    estimate=(w-np.eye(2))/(h*h)
    # Our Wilson convention dW=-A.dR W gives W_loop = I - F_xy area + O(h^3).
    assert np.linalg.norm(estimate+f)/max(np.linalg.norm(f),1e-30) < 1.5e-2


def test_reverse_path_is_inverse_transport():
    field=variable_field(); p=dogleg_path_v293(np.array([-.2,.1]),np.array([.5,-.3]),axis=1)
    w=wilson_path_v293(field,p); wr=wilson_path_v293(field,p.reversed())
    assert np.linalg.norm(wr@w-np.eye(2)) < 2e-10


def test_patch_transition_matches_single_patch_transport():
    source=variable_field(); gauge=coordinate_gauge(); target=gauge_transform_connection_v293(source,gauge)
    transition=PatchTransitionV293("alpha","beta",gauge)
    probes=[np.array([-.2,.1]),np.array([.1,-.15]),np.array([.4,.2])]
    assert validate_patch_transition_v293(source,target,transition,probes) < 3e-7
    a=np.array([-.45,-.10]); x=np.array([.02,.18]); b=np.array([.42,.35])
    mixed=two_patch_transport_v293(source,target,transition,a,x,b)
    same_path=PolylinePathV293(np.vstack([a,x,b]),label="same-two-segment-path").validate()
    reference=gauge.at(b).conj().T@wilson_path_v293(source,same_path)
    assert np.linalg.norm(mixed-reference) < 2e-9


def synthetic_snapshot(state,field,quad,candidate):
    j=finite_difference_full_jacobian_v293(state,field,quad,relative_step=4e-6)
    k=candidate_coefficient_block_v293(candidate,state.nstate,field,quad)
    residual=0.13*k[:,0]-0.09*k[:,1]+0.025j*k[:,0]
    jr=np.vstack([j.real,j.imag]); rr=np.concatenate([residual.real,residual.imag])
    correction=np.linalg.lstsq(jr,-rr,rcond=1e-10)[0]; rr=rr+jr@correction
    n=j.shape[0]; residual=rr[:n]+1j*rr[n:]
    m=np.real(j.conj().T@j); s=np.linalg.svd(m,compute_uv=False)
    rec=SVDReceiptV292(s,1e-12,len(s),float(s[0]/s[-1]),0.0,0.0)
    return CurvedTDVPSnapshotV292(j,m,np.zeros(j.shape[1]),np.zeros(j.shape[1]),residual,rec).validate()


def test_v292_residual_spawning_rule_survives_coordinate_dependent_connection():
    state=fixture_state(); field=variable_field(); quad=tensor_legendre_quadrature_v292(((-2.0,2.0),(-2.0,2.0)),(4,4))
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    snap=synthetic_snapshot(state,field,quad,candidate)
    settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
    ev=evaluate_curved_spawn_candidate_v293(state,field,quad,snap,candidate,settings=settings)
    k=candidate_coefficient_block_v293(candidate,state.nstate,field,quad)
    actual=actual_linear_residual_reduction_v292(snap,k)
    assert ev.admitted and actual>0
    assert abs(ev.predicted_residual_reduction-actual)/actual < 3e-6


def test_spawn_projection_remains_wavefunction_preserving_with_variable_connection():
    state=fixture_state(); field=variable_field(); quad=tensor_legendre_quadrature_v292(((-2.0,2.0),(-2.0,2.0)),(5,5))
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    out,receipt=spawn_curved_candidate_v293(state,candidate,field,quad)
    assert out.ngaussian==2
    assert receipt.relative_loss < 1e-20
    assert receipt.fidelity > 1-1e-12
    assert np.linalg.norm(out.coefficients[-1]) < 1e-9


def test_packet_section_is_pointwise_gauge_covariant_when_center_coefficient_is_transformed():
    state=fixture_state(); _,source=base_constant(); gauge=coordinate_gauge(); target=gauge_transform_connection_v293(source,gauge)
    quad=tensor_legendre_quadrature_v292(((-1.6,1.6),(-1.6,1.6)),(4,4))
    psi_source=np.einsum('nsia,ia->ns',coefficient_basis_v293(state,source,quad),state.coefficients,optimize=True)
    c2=state.coefficients.copy(); c2[0]=gauge.at(state.q[0]).conj().T@c2[0]
    transformed_state=CorrelatedGaussianSpinorStateV270(state.q,state.p,state.width_matrices,state.chirp_matrices,c2,state.time_au).validate(require_normalized=False)
    psi_target=np.einsum('nsia,ia->ns',coefficient_basis_v293(transformed_state,target,quad),transformed_state.coefficients,optimize=True)
    expected=np.asarray([gauge.at(r).conj().T@psi_source[n] for n,r in enumerate(quad.points)])
    assert np.max(np.abs(psi_target-expected)) < 4e-12


def test_full_spawn_event_runs_under_coordinate_dependent_connection():
    from gaussian_dynamics.curved_basis_adaptation_v293 import adapt_curved_basis_once_v293
    state=fixture_state(); field=variable_field(); quad=tensor_legendre_quadrature_v292(((-1.8,1.8),(-1.8,1.8)),(4,4))
    candidate=generate_multidimensional_spawn_candidates_v270(state)[0]
    snap=synthetic_snapshot(state,field,quad,candidate)
    settings=CurvedBasisSettingsV292(minimum_predicted_residual_reduction=1e-14,minimum_relative_residual_reduction=1e-12)
    def identity_h(state,psi,quad): return np.asarray(psi,complex)
    event=adapt_curved_basis_once_v293(state,field,quad,snap,identity_h,settings=settings)
    assert event.event_kind=='spawn'
    assert event.after.ngaussian==2
    assert event.active_shape_mask_after.tolist()==[True,False]
    assert event.projection.relative_loss < 1e-18
