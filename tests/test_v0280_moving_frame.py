from dataclasses import replace
import numpy as np
import pytest
from gaussian_dynamics.correlated_basis_adaptation_v270 import adapt_multidimensional_basis_once_v270
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270,_velocity_blocks_v270,build_correlated_metric_system_v270,correlated_variational_energy_v270,evaluate_correlated_state_v270
from gaussian_dynamics.moving_frame_v280 import FlatMovingFrameV280,MovingFrameCorrelatedStateV280,adapt_moving_frame_basis_once_v280,evaluate_moving_physical_v280,evaluate_moving_section_v280,fixed_to_moving_state_v280,moving_frame_hamiltonian_v280,moving_frame_implicit_midpoint_step_v280,moving_frame_velocity_v280,moving_to_fixed_state_v280,reference_wavefunction_error_v280,require_flat_moving_frame_v280
from gaussian_dynamics.multidimensional_soc_v260 import two_state_ci_soc_model_v260

def _frame():
    return FlatMovingFrameV280(np.asarray([[0.23,0.31-0.17j],[0.31+0.17j,-0.11]],complex),np.asarray([0.37,-0.22]),np.asarray([[0.08,-0.03],[-0.03,-0.05]]),0.19).validate()
def _fixed_state():
    return CorrelatedGaussianSpinorStateV270([[-0.25,0.1]],[[3.0,0.2]],[[[1.7,0.25],[0.25,2.6]]],[[[0.02,0.03],[0.03,-0.01]]],[[0.92+0.04j,0.18-0.09j]]).normalized()
def _moving_state(): return fixed_to_moving_state_v280(_fixed_state(),_frame()).validate(_frame(),require_normalized=True)
def _constant_unitary():
    angle=0.431; phase=np.exp(0.37j); return np.asarray([[np.cos(angle),-phase.conjugate()*np.sin(angle)],[phase*np.sin(angle),np.cos(angle)]],complex)

def test_v0280_frame_unitary_at_arbitrary_points():
    gauges=_frame().unitary(np.asarray([[-0.7,0.2],[0.1,-0.4],[0.8,0.5]])); assert max(np.linalg.norm(g.conj().T@g-np.eye(2)) for g in gauges)<3e-15
def test_v0280_connection_is_antihermitian():
    c=_frame().connection([0.17,-0.31]); assert np.max(np.abs(c+c.conj().swapaxes(-1,-2)))<3e-15
def test_v0280_exact_curvature_is_zero(): assert np.max(np.abs(_frame().curvature([0.22,-0.14])))==0.0
def test_v0280_transporter_is_unitary():
    W=_frame().transporter([0.43,-0.17],[-0.21,0.33]); assert np.linalg.norm(W.conj().T@W-np.eye(2))<4e-15
def test_v0280_transporter_is_identity_at_center():
    q=np.asarray([0.31,-0.27]); assert np.max(np.abs(_frame().transporter(q,q)-np.eye(2)))<3e-15
def test_v0280_transporter_composes_exactly():
    f=_frame(); a=np.asarray([-0.3,0.1]); b=np.asarray([0.2,-0.4]); c=np.asarray([0.7,0.25]); assert np.max(np.abs(f.transporter(c,a)-f.transporter(c,b)@f.transporter(b,a)))<1.5e-15
def test_v0280_reference_moving_state_roundtrip():
    fixed=_fixed_state(); rec=moving_to_fixed_state_v280(fixed_to_moving_state_v280(fixed,_frame()),_frame()); assert np.max(np.abs(rec.coefficients-fixed.coefficients))<6e-16; assert np.array_equal(rec.q,fixed.q); assert np.array_equal(rec.width_matrices,fixed.width_matrices)
def test_v0280_physical_wavefunction_matches_fixed_reference():
    points=np.asarray([[-0.5,-0.1],[-0.2,0.15],[0.4,-0.3],[0.9,0.2]]); assert reference_wavefunction_error_v280(_moving_state(),_frame(),points)<8e-16
def test_v0280_explicit_physical_evaluation_matches_fixed_evaluator():
    f=_frame(); s=_moving_state(); p=np.asarray([[0.11,-0.27],[0.62,0.09]]); assert np.max(np.abs(evaluate_moving_physical_v280(s,f,p)-evaluate_correlated_state_v270(moving_to_fixed_state_v280(s,f),p)))<8e-16
def test_v0280_local_section_transforms_under_constant_gauge():
    f=_frame();s=_moving_state();U=_constant_unitary();p=np.asarray([[-0.2,0.3],[0.5,-0.1]]); base=evaluate_moving_section_v280(s,f,p); trans=evaluate_moving_section_v280(s.gauge_rotated(U),f.gauge_rotated(U),p); assert np.max(np.abs(trans-np.einsum('ab,...b->...a',U.conj().T,base)))<9e-16
def test_v0280_physical_wavefunction_is_constant_gauge_invariant():
    f=_frame();s=_moving_state();U=_constant_unitary();p=np.asarray([[-0.2,0.3],[0.5,-0.1]]); assert np.max(np.abs(evaluate_moving_physical_v280(s,f,p)-evaluate_moving_physical_v280(s.gauge_rotated(U),f.gauge_rotated(U),p)))<9e-16
def test_v0280_transporter_is_constant_gauge_covariant():
    f=_frame();U=_constant_unitary();p=np.asarray([0.41,-0.22]);q=np.asarray([-0.19,0.16]);base=f.transporter(p,q); assert np.max(np.abs(f.gauge_rotated(U).transporter(p,q)-U.conj().T@base@U))<8e-16
def test_v0280_connection_is_constant_gauge_covariant():
    f=_frame();U=_constant_unitary();p=np.asarray([0.41,-0.22]);base=f.connection(p); expected=np.einsum('ab,xbc,cd->xad',U.conj().T,base,U); assert np.max(np.abs(f.gauge_rotated(U).connection(p)-expected))<9e-16
def test_v0280_moving_hamiltonian_is_hermitian():
    value=moving_frame_hamiltonian_v280(two_state_ci_soc_model_v260(),_frame(),[[-0.2,0.3],[0.5,-0.1]]); assert np.max(np.abs(value-value.conj().swapaxes(-1,-2)))<3e-15
def test_v0280_moving_hamiltonian_is_constant_gauge_covariant():
    m=two_state_ci_soc_model_v260();f=_frame();U=_constant_unitary();p=np.asarray([0.12,-0.33]);base=moving_frame_hamiltonian_v280(m,f,p); assert np.max(np.abs(moving_frame_hamiltonian_v280(m,f.gauge_rotated(U),p)-U.conj().T@base@U))<9e-16
def test_v0280_velocity_maps_back_to_fixed_tdvp_tangent():
    f=_frame();mov=_moving_state();fixed=moving_to_fixed_state_v280(mov,f); vel,sys=moving_frame_velocity_v280(mov,two_state_ci_soc_model_v260(),f); cmdot,qdot,_,_,_=_velocity_blocks_v270(fixed,vel); eps=2e-7; pert=MovingFrameCorrelatedStateV280(mov.q+eps*qdot,mov.p,mov.width_matrices,mov.chirp_matrices,mov.center_coefficients+eps*cmdot,mov.time_au); finite=(moving_to_fixed_state_v280(pert,f).coefficients-fixed.coefficients)/eps; fixed_cdot=_velocity_blocks_v270(fixed,sys.velocity)[0]; assert np.max(np.abs(finite-fixed_cdot))<2e-7
def test_v0280_velocity_physical_reference_matches_v0270():
    f=_frame();mov=_moving_state();fixed=moving_to_fixed_state_v280(mov,f);_,system=moving_frame_velocity_v280(mov,two_state_ci_soc_model_v260(),f);direct=build_correlated_metric_system_v270(fixed,two_state_ci_soc_model_v260());assert np.max(np.abs(system.velocity-direct.velocity))<1e-13
def test_v0280_midpoint_endpoint_maps_exactly_to_v0270():
    f=_frame();step=moving_frame_implicit_midpoint_step_v280(_moving_state(),two_state_ci_soc_model_v260(),f,0.001);fixed_end=moving_to_fixed_state_v280(step.end,f);assert np.max(np.abs(fixed_end.coefficients-step.fixed_step.end.coefficients))<8e-16;assert np.max(np.abs(fixed_end.q-step.fixed_step.end.q))==0.0
def test_v0280_midpoint_preserves_v0270_norm_and_energy_receipts():
    f=_frame();s=_moving_state();m=two_state_ci_soc_model_v260();step=moving_frame_implicit_midpoint_step_v280(s,m,f,0.001);assert abs(step.fixed_step.norm_change)<1e-8;assert abs(step.fixed_step.energy_change_hartree)<1e-8;assert abs(correlated_variational_energy_v270(moving_to_fixed_state_v280(step.end,f),m)-step.fixed_step.end_energy_hartree)<2e-13
def test_v0280_midpoint_physical_endpoint_is_constant_gauge_invariant():
    f=_frame();s=_moving_state();U=_constant_unitary();m=two_state_ci_soc_model_v260();base=moving_frame_implicit_midpoint_step_v280(s,m,f,0.001);trans=moving_frame_implicit_midpoint_step_v280(s.gauge_rotated(U),m,f.gauge_rotated(U),0.001);fb=moving_to_fixed_state_v280(base.end,f);ft=moving_to_fixed_state_v280(trans.end,f.gauge_rotated(U));assert np.max(np.abs(fb.coefficients-ft.coefficients))<2e-13;assert np.max(np.abs(fb.q-ft.q))<2e-13
def test_v0280_lifecycle_event_matches_fixed_trivialization():
    f=_frame();mov=_moving_state();m=two_state_ci_soc_model_v260();fixed_event=adapt_multidimensional_basis_once_v270(moving_to_fixed_state_v280(mov,f),m);moving_event=adapt_moving_frame_basis_once_v280(mov,m,f);assert moving_event.event_type==fixed_event.event_kind;assert moving_event.reason==fixed_event.reason;assert np.max(np.abs(moving_to_fixed_state_v280(moving_event.after,f).coefficients-fixed_event.after.coefficients))<8e-16
def test_v0280_nonunitary_right_gauge_is_rejected():
    with pytest.raises(ValueError,match='unitary'): replace(_frame(),right_unitary=np.asarray([[1.0,0.2],[0,1.0]])).validate()
def test_v0280_nonhermitian_generator_is_rejected():
    f=_frame();broken=f.generator.copy();broken[0,1]+=0.2
    with pytest.raises(ValueError,match='Hermitian'): replace(f,generator=broken).validate()
def test_v0280_nonsymmetric_phase_hessian_is_rejected():
    with pytest.raises(ValueError,match='symmetric'): replace(_frame(),phase_hessian=np.asarray([[0.1,0.2],[0.0,0.1]])).validate()
def test_v0280_nonflat_connection_is_rejected():
    curvature=np.zeros((2,2,2,2),complex);curvature[0,1,0,1]=1e-3
    with pytest.raises(ValueError,match='non-flat'): require_flat_moving_frame_v280(_frame(),curvature=curvature)
def test_v0280_missing_trivialization_is_rejected():
    with pytest.raises(ValueError,match='trivialization'): require_flat_moving_frame_v280(_frame(),trivialization_available=False)
def test_v0280_state_frame_dimension_mismatch_is_rejected():
    wrong=FlatMovingFrameV280(np.eye(2),[0.2],[[0.0]],0.0).validate()
    with pytest.raises(ValueError,match='dimension'): _moving_state().validate(wrong)
def test_v0280_nonunitary_state_gauge_rotation_is_rejected():
    with pytest.raises(ValueError,match='unitary'): _moving_state().gauge_rotated([[1.0,0.2],[0,1.0]])
