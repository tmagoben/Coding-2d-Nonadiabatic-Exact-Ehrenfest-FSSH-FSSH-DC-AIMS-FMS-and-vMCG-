import numpy as np
import pytest
from gaussian_dynamics.moving_frame_v280 import FlatMovingFrameV280
from gaussian_dynamics.moving_frame_validation_v280 import build_lattice_gauge_oracle_v280,finite_difference_connection_residual_v280,finite_difference_curvature_residual_v280,lattice_action_covariance_v280,lattice_propagation_covariance_v280,periodic_second_derivative_v280,scalar_lattice_kinetic_v280
from gaussian_dynamics.multidimensional_soc_v260 import QuadraticSpinHamiltonianNDV260,two_state_ci_soc_model_v260

def _frame(): return FlatMovingFrameV280(np.asarray([[0.2,0.3-0.1j],[0.3+0.1j,-0.15]]),[0.31,-0.27],[[0.06,0.02],[0.02,-0.04]],0.13).validate()
def _axes(): return (np.linspace(-0.8,0.8,4,endpoint=False),np.linspace(-0.6,0.6,3,endpoint=False))
def _oracle(): return build_lattice_gauge_oracle_v280(two_state_ci_soc_model_v260(),_frame(),_axes())
def test_v0280_connection_matches_centered_finite_difference(): assert finite_difference_connection_residual_v280(_frame(),[0.17,-0.29])<2e-10
def test_v0280_curvature_matches_zero_by_independent_finite_difference(): assert finite_difference_curvature_residual_v280(_frame(),[0.17,-0.29])<3e-11
def test_v0280_periodic_second_derivative_is_symmetric_and_annihilates_constant():
    d2=periodic_second_derivative_v280(5,0.2); assert np.max(np.abs(d2-d2.T))==0.0; assert np.max(np.abs(d2@np.ones(5)))<2e-14
def test_v0280_lattice_hamiltonians_are_hermitian(): assert _oracle().hamiltonian_hermiticity_residual<3e-15
def test_v0280_lattice_link_transport_is_unitary(): assert _oracle().maximum_link_unitarity_residual<3e-15
def test_v0280_lattice_direct_links_equal_global_unitary_similarity(): assert _oracle().similarity_residual<8e-16
def test_v0280_lattice_action_is_gauge_covariant():
    o=_oracle();rng=np.random.default_rng(280);s=rng.normal(size=o.fixed_hamiltonian.shape[0])+1j*rng.normal(size=o.fixed_hamiltonian.shape[0]);s/=np.linalg.norm(s);assert lattice_action_covariance_v280(o,s)<8e-16
def test_v0280_lattice_propagation_is_gauge_covariant():
    o=_oracle();rng=np.random.default_rng(281);s=rng.normal(size=o.fixed_hamiltonian.shape[0])+1j*rng.normal(size=o.fixed_hamiltonian.shape[0]);s/=np.linalg.norm(s);assert lattice_propagation_covariance_v280(o,s,0.07)<2e-15
def test_v0280_lattice_transformation_is_unitary():
    o=_oracle();assert np.linalg.norm(o.transformation.conj().T@o.transformation-np.eye(o.transformation.shape[0]))<5e-15
def test_v0280_lattice_oracle_rejects_nondiagonal_mass_matrix():
    base=two_state_ci_soc_model_v260();mass=base.mass_matrix_au.copy();mass[0,1]=mass[1,0]=3.0;model=QuadraticSpinHamiltonianNDV260(mass,base.H0,base.H1,base.H2,label=base.label,model_kind=base.model_kind,projectors=base.projectors,complete_spin_manifold=base.complete_spin_manifold,physical_soc=base.physical_soc,soc_scale_hartree=base.soc_scale_hartree,source=base.source).validate()
    with pytest.raises(ValueError,match='diagonal mass'): scalar_lattice_kinetic_v280(_axes(),model.mass_matrix_au)
