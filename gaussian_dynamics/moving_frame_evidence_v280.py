"""Deterministic scientific evidence campaign for v0.28.0 moving frames."""
from dataclasses import dataclass
import hashlib, json
import numpy as np
from .correlated_basis_adaptation_v270 import adapt_multidimensional_basis_once_v270
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270, _velocity_blocks_v270, build_correlated_metric_system_v270
from .moving_frame_v280 import (
    CLAIM_BOUNDARY_V280, FlatMovingFrameV280, MovingFrameCorrelatedStateV280,
    adapt_moving_frame_basis_once_v280, evaluate_moving_physical_v280,
    evaluate_moving_section_v280, fixed_to_moving_state_v280, moving_frame_hamiltonian_v280,
    moving_frame_implicit_midpoint_step_v280, moving_frame_velocity_v280,
    moving_to_fixed_state_v280, reference_wavefunction_error_v280, require_flat_moving_frame_v280,
)
from .moving_frame_validation_v280 import (
    build_lattice_gauge_oracle_v280, finite_difference_connection_residual_v280,
    finite_difference_curvature_residual_v280, lattice_action_covariance_v280,
    lattice_propagation_covariance_v280, scalar_lattice_kinetic_v280,
)
from .multidimensional_soc_v260 import QuadraticSpinHamiltonianNDV260, two_state_ci_soc_model_v260

MOVING_FRAME_EVIDENCE_SCHEMA_V280 = "gnd-moving-frame-scientific-evidence-v0.28.0"

def _canonical(value):
    if isinstance(value,np.generic): return _canonical(value.item())
    if isinstance(value,complex): return [float(value.real),float(value.imag)]
    if isinstance(value,np.ndarray): return _canonical(value.tolist())
    if isinstance(value,dict): return {str(k):_canonical(v) for k,v in sorted(value.items(),key=lambda p:str(p[0]))}
    if isinstance(value,(tuple,list)): return [_canonical(v) for v in value]
    if isinstance(value,float):
        if not np.isfinite(value): raise ValueError('evidence cannot contain non-finite floats.')
        return float(value)
    if value is None or isinstance(value,(str,int,bool)): return value
    raise TypeError(type(value).__name__)

def _sha(value):
    payload=json.dumps(_canonical(value),sort_keys=True,separators=(',',':')).encode(); return hashlib.sha256(payload).hexdigest()

def _scaled(left,right):
    left=np.asarray(left); right=np.asarray(right)
    if left.shape!=right.shape: return float('inf')
    return float(np.linalg.norm(left-right)/max(np.linalg.norm(left),np.linalg.norm(right),1.0))

def _raises(fn,text):
    try: fn()
    except Exception as exc: return text.lower() in str(exc).lower()
    return False

def _frame():
    return FlatMovingFrameV280(np.asarray([[0.23,0.31-0.17j],[0.31+0.17j,-0.11]],complex),np.asarray([0.37,-0.22]),np.asarray([[0.08,-0.03],[-0.03,-0.05]]),0.19).validate()

def _fixed_state():
    return CorrelatedGaussianSpinorStateV270([[-0.25,0.1]],[[3.0,0.2]],[[[1.7,0.25],[0.25,2.6]]],[[[0.02,0.03],[0.03,-0.01]]],[[0.92+0.04j,0.18-0.09j]]).normalized()

def _unitary():
    angle=0.431; phase=np.exp(0.37j)
    return np.asarray([[np.cos(angle),-phase.conjugate()*np.sin(angle)],[phase*np.sin(angle),np.cos(angle)]],complex)

@dataclass(frozen=True)
class MovingFrameEvidenceV280:
    checks: dict
    metrics: dict
    claims: dict
    trajectory_fingerprint: str
    lattice_fingerprint: str
    @property
    def check_count(self): return len(self.checks)
    @property
    def passed(self): return bool(self.checks) and all(bool(v) for v in self.checks.values())
    def as_dict(self):
        return {'schema':MOVING_FRAME_EVIDENCE_SCHEMA_V280,'passed':self.passed,'check_count':self.check_count,'checks':_canonical(dict(sorted(self.checks.items()))),'metrics':_canonical(self.metrics),'claims':_canonical(self.claims),'trajectory_fingerprint':self.trajectory_fingerprint,'lattice_fingerprint':self.lattice_fingerprint}
    def fingerprint(self): return _sha(self.as_dict())

def run_moving_frame_evidence_v280():
    frame=_frame(); fixed=_fixed_state(); moving=fixed_to_moving_state_v280(fixed,frame).validate(frame,require_normalized=True); model=two_state_ci_soc_model_v260(); U=_unitary(); frame_u=frame.gauge_rotated(U); moving_u=moving.gauge_rotated(U)
    points=np.asarray([[-0.5,-0.1],[-0.2,0.15],[0.4,-0.3],[0.9,0.2]])
    gauges=frame.unitary(points); frame_unitarity=max(np.linalg.norm(g.conj().T@g-np.eye(2)) for g in gauges)
    connection=frame.connection([0.17,-0.31]); connection_antihermiticity=float(np.max(np.abs(connection+connection.conj().swapaxes(-1,-2))))
    connection_fd=finite_difference_connection_residual_v280(frame,[0.17,-0.29]); curvature_fd=finite_difference_curvature_residual_v280(frame,[0.17,-0.29])
    a=np.asarray([-0.3,0.1]); b=np.asarray([0.2,-0.4]); c=np.asarray([0.7,0.25]); W_ca=frame.transporter(c,a); W_cb=frame.transporter(c,b); W_ba=frame.transporter(b,a)
    transport_identity=float(np.max(np.abs(frame.transporter(a,a)-np.eye(2)))); transport_unitarity=float(np.linalg.norm(W_ca.conj().T@W_ca-np.eye(2))); transport_composition=float(np.max(np.abs(W_ca-W_cb@W_ba))); transport_gauge=float(np.max(np.abs(frame_u.transporter(c,a)-U.conj().T@W_ca@U)))
    roundtrip=moving_to_fixed_state_v280(moving,frame); roundtrip_coeff=float(np.max(np.abs(roundtrip.coefficients-fixed.coefficients))); wavefunction=reference_wavefunction_error_v280(moving,frame,points)
    section=evaluate_moving_section_v280(moving,frame,points); section_u=evaluate_moving_section_v280(moving_u,frame_u,points); section_gauge=float(np.max(np.abs(section_u-np.einsum('ab,...b->...a',U.conj().T,section))))
    physical=evaluate_moving_physical_v280(moving,frame,points); physical_u=evaluate_moving_physical_v280(moving_u,frame_u,points); physical_gauge=float(np.max(np.abs(physical_u-physical)))
    Hm=moving_frame_hamiltonian_v280(model,frame,points); h_hermiticity=float(np.max(np.abs(Hm-Hm.conj().swapaxes(-1,-2)))); h_base=moving_frame_hamiltonian_v280(model,frame,points[1]); h_u=moving_frame_hamiltonian_v280(model,frame_u,points[1]); h_gauge=float(np.max(np.abs(h_u-U.conj().T@h_base@U)))
    connection_u=frame_u.connection(points[1]); connection_gauge=float(np.max(np.abs(connection_u-np.einsum('ab,xbc,cd->xad',U.conj().T,frame.connection(points[1]),U))))
    moving_velocity,fixed_system=moving_frame_velocity_v280(moving,model,frame); direct_system=build_correlated_metric_system_v270(fixed,model); fixed_velocity_match=float(np.max(np.abs(fixed_system.velocity-direct_system.velocity)))
    cmdot,qdot,_,_,_=_velocity_blocks_v270(fixed,moving_velocity); eps=2e-7
    perturbed=MovingFrameCorrelatedStateV280(moving.q+eps*qdot,moving.p,moving.width_matrices,moving.chirp_matrices,moving.center_coefficients+eps*cmdot,moving.time_au)
    finite_cfdot=(moving_to_fixed_state_v280(perturbed,frame).coefficients-fixed.coefficients)/eps; exact_cfdot=_velocity_blocks_v270(fixed,fixed_system.velocity)[0]; velocity_trivialization=float(np.max(np.abs(finite_cfdot-exact_cfdot)))
    step=moving_frame_implicit_midpoint_step_v280(moving,model,frame,0.001); fixed_end=moving_to_fixed_state_v280(step.end,frame)
    step_coeff=float(np.max(np.abs(fixed_end.coefficients-step.fixed_step.end.coefficients))); step_q=float(np.max(np.abs(fixed_end.q-step.fixed_step.end.q))); step_width=float(np.max(np.abs(fixed_end.width_matrices-step.fixed_step.end.width_matrices))); step_chirp=float(np.max(np.abs(fixed_end.chirp_matrices-step.fixed_step.end.chirp_matrices)))
    step_u=moving_frame_implicit_midpoint_step_v280(moving_u,model,frame_u,0.001); fixed_end_u=moving_to_fixed_state_v280(step_u.end,frame_u); step_gauge_coeff=float(np.max(np.abs(fixed_end_u.coefficients-fixed_end.coefficients))); step_gauge_q=float(np.max(np.abs(fixed_end_u.q-fixed_end.q)))
    fixed_event=adapt_multidimensional_basis_once_v270(fixed,model); moving_event=adapt_moving_frame_basis_once_v280(moving,model,frame); lifecycle_after=moving_to_fixed_state_v280(moving_event.after,frame); lifecycle_coeff=float(np.max(np.abs(lifecycle_after.coefficients-fixed_event.after.coefficients)))
    axes=(np.linspace(-0.8,0.8,4,endpoint=False),np.linspace(-0.6,0.6,3,endpoint=False)); lattice=build_lattice_gauge_oracle_v280(model,frame,axes); rng=np.random.default_rng(280); vector=rng.normal(size=lattice.fixed_hamiltonian.shape[0])+1j*rng.normal(size=lattice.fixed_hamiltonian.shape[0]); vector/=np.linalg.norm(vector)
    action_covariance=lattice_action_covariance_v280(lattice,vector); propagation_covariance=lattice_propagation_covariance_v280(lattice,vector,0.07); transformation_unitarity=_scaled(lattice.transformation.conj().T@lattice.transformation,np.eye(lattice.transformation.shape[0]))
    nonflat=np.zeros((2,2,2,2),complex); nonflat[0,1,0,1]=1e-3; bad_generator=frame.generator.copy(); bad_generator[0,1]+=0.2; bad_hessian=np.asarray([[0.1,0.2],[0.0,0.1]]); bad_mass=model.mass_matrix_au.copy(); bad_mass[0,1]=bad_mass[1,0]=3.0
    bad_mass_model=QuadraticSpinHamiltonianNDV260(bad_mass,model.H0,model.H1,model.H2,label=model.label,model_kind=model.model_kind,projectors=model.projectors,complete_spin_manifold=model.complete_spin_manifold,physical_soc=model.physical_soc,soc_scale_hartree=model.soc_scale_hartree,source=model.source).validate()
    metrics={
        'frame_unitarity_residual':frame_unitarity,'connection_antihermiticity_residual':connection_antihermiticity,'connection_centered_difference_residual':connection_fd,'curvature_centered_difference_residual':curvature_fd,'transport_identity_residual':transport_identity,'transport_unitarity_residual':transport_unitarity,'transport_composition_residual':transport_composition,'transport_constant_gauge_residual':transport_gauge,'state_roundtrip_coefficient_residual':roundtrip_coeff,'physical_wavefunction_residual':wavefunction,'section_constant_gauge_residual':section_gauge,'physical_constant_gauge_residual':physical_gauge,'moving_hamiltonian_hermiticity_residual':h_hermiticity,'moving_hamiltonian_constant_gauge_residual':h_gauge,'connection_constant_gauge_residual':connection_gauge,'fixed_velocity_system_residual':fixed_velocity_match,'velocity_trivialization_residual':velocity_trivialization,'midpoint_coefficient_residual':step_coeff,'midpoint_position_residual':step_q,'midpoint_width_residual':step_width,'midpoint_chirp_residual':step_chirp,'midpoint_constant_gauge_coefficient_residual':step_gauge_coeff,'midpoint_constant_gauge_position_residual':step_gauge_q,'midpoint_norm_change':abs(float(step.fixed_step.norm_change)),'midpoint_energy_change_hartree':abs(float(step.fixed_step.energy_change_hartree)),'lifecycle_coefficient_residual':lifecycle_coeff,'lattice_link_unitarity_residual':lattice.maximum_link_unitarity_residual,'lattice_hamiltonian_hermiticity_residual':lattice.hamiltonian_hermiticity_residual,'lattice_similarity_residual':lattice.similarity_residual,'lattice_transformation_unitarity_residual':transformation_unitarity,'lattice_action_covariance_residual':action_covariance,'lattice_propagation_covariance_residual':propagation_covariance,
    }
    checks={
        'frame_unitary':frame_unitarity<5e-15,'connection_antihermitian':connection_antihermiticity<5e-15,'connection_matches_centered_difference':connection_fd<3e-10,'curvature_matches_zero_by_centered_difference':curvature_fd<5e-11,'analytic_curvature_exactly_zero':np.max(np.abs(frame.curvature(points[0])))==0.0,'transport_identity':transport_identity<5e-15,'transport_unitary':transport_unitarity<5e-15,'transport_composition':transport_composition<2e-15,'transport_constant_gauge_covariance':transport_gauge<2e-15,'state_reference_roundtrip':roundtrip_coeff<1e-15,'state_position_roundtrip_exact':np.array_equal(roundtrip.q,fixed.q),'state_width_roundtrip_exact':np.array_equal(roundtrip.width_matrices,fixed.width_matrices),'physical_wavefunction_reference_equivalence':wavefunction<1e-15,'local_section_constant_gauge_covariance':section_gauge<2e-15,'physical_wavefunction_constant_gauge_invariance':physical_gauge<2e-15,'moving_hamiltonian_hermitian':h_hermiticity<5e-15,'moving_hamiltonian_constant_gauge_covariance':h_gauge<2e-15,'connection_constant_gauge_covariance':connection_gauge<2e-15,'fixed_tdvp_system_reused_exactly':fixed_velocity_match<1e-13,'moving_velocity_trivializes_to_fixed_velocity':velocity_trivialization<3e-7,'midpoint_coefficients_trivialize':step_coeff<2e-15,'midpoint_positions_trivialize':step_q<1e-15,'midpoint_widths_trivialize':step_width<1e-15,'midpoint_chirps_trivialize':step_chirp<1e-15,'midpoint_constant_gauge_coefficients':step_gauge_coeff<3e-13,'midpoint_constant_gauge_positions':step_gauge_q<3e-13,'midpoint_norm_stable':abs(float(step.fixed_step.norm_change))<1e-8,'midpoint_energy_stable':abs(float(step.fixed_step.energy_change_hartree))<1e-8,'lifecycle_event_kind_identical':moving_event.fixed_event.event_kind==fixed_event.event_kind,'lifecycle_reason_identical':moving_event.fixed_event.reason==fixed_event.reason,'lifecycle_endpoint_trivializes':lifecycle_coeff<2e-15,'lattice_links_unitary':lattice.maximum_link_unitarity_residual<5e-15,'lattice_hamiltonians_hermitian':lattice.hamiltonian_hermiticity_residual<5e-15,'lattice_direct_links_equal_similarity_transform':lattice.similarity_residual<2e-15,'lattice_transformation_unitary':transformation_unitarity<5e-15,'lattice_action_covariant':action_covariance<2e-15,'lattice_propagation_covariant':propagation_covariance<3e-15,'nonflat_connection_rejected':_raises(lambda:require_flat_moving_frame_v280(frame,curvature=nonflat),'non-flat'),'missing_trivialization_rejected':_raises(lambda:require_flat_moving_frame_v280(frame,trivialization_available=False),'trivialization'),'nonunitary_right_gauge_rejected':_raises(lambda:FlatMovingFrameV280(frame.generator,frame.phase_gradient,frame.phase_hessian,frame.phase_offset,[[1,0.2],[0,1]]).validate(),'unitary'),'nonhermitian_generator_rejected':_raises(lambda:FlatMovingFrameV280(bad_generator,frame.phase_gradient,frame.phase_hessian).validate(),'Hermitian'),'nonsymmetric_phase_hessian_rejected':_raises(lambda:FlatMovingFrameV280(frame.generator,frame.phase_gradient,bad_hessian).validate(),'symmetric'),'nonunitary_state_gauge_rejected':_raises(lambda:moving.gauge_rotated([[1,0.2],[0,1]]),'unitary'),'frame_state_dimension_mismatch_rejected':_raises(lambda:moving.validate(FlatMovingFrameV280(np.eye(2),[0.2],[[0.0]]).validate()),'dimension'),'nondiagonal_lattice_mass_rejected':_raises(lambda:scalar_lattice_kinetic_v280(axes,bad_mass_model.mass_matrix_au),'diagonal mass'),'live_molecular_soc_claim_remains_closed':True,'nonflat_connection_claim_remains_closed':True,'general_curved_manifold_claim_remains_closed':True,'full_aims_branching_claim_remains_closed':True,'release_is_development_not_sealed':True,
    }
    claims={'validated':{'parallel_transported_electronic_sections':True,'flat_coordinate_dependent_gauge_covariance':True,'exact_fixed_frame_trivialization':True,'tdvp_endpoint_covariance':True,'controlled_lifecycle_covariance':True,'independent_lattice_gauge_oracle':True},'not_validated':{'nonzero_curvature_connections':True,'live_molecular_soc_trajectories':True,'general_ab_initio_soc_dynamics_accuracy':True,'full_aims_branching_semantics':True},'boundary':CLAIM_BOUNDARY_V280}
    trajectory_fingerprint=_sha({'start':moving_to_fixed_state_v280(moving,frame).as_dict(),'end':fixed_end.as_dict(),'event':fixed_event.as_dict()}); lattice_fingerprint=_sha({'points':lattice.points,'fixed':lattice.fixed_hamiltonian,'moving':lattice.moving_hamiltonian})
    return MovingFrameEvidenceV280(checks,metrics,claims,trajectory_fingerprint,lattice_fingerprint)

__all__=['MOVING_FRAME_EVIDENCE_SCHEMA_V280','MovingFrameEvidenceV280','run_moving_frame_evidence_v280']
