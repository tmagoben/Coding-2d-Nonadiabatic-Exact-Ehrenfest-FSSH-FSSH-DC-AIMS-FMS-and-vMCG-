import numpy as np

from gaussian_dynamics.grids import uniform_grid
from gaussian_dynamics.gaussian import frozen_gaussian
from gaussian_dynamics.potentials import avoided_crossing_diabatic
from gaussian_dynamics.exact_multistate import run_multistate_exact
from gaussian_dynamics.adiabatic import (
    adiabatic_grid, adiabatic_hamiltonian_action,
    diabatic_hamiltonian_action, transform_ad_to_dia, transform_dia_to_ad
)
from gaussian_dynamics.adiabatic_spawning import (
    AdiabaticTBF, spawn_child_energy_conserving, maybe_spawn,
    adiabatic_gaussian_basis_matrices, expand_coefficients_for_new_basis
)


def test_multistate_exact_norm():
    x,dx=uniform_grid(-10,10,512)
    V=avoided_crossing_diabatic(x)
    psi=np.zeros((len(x),2),complex)
    psi[:,0]=frozen_gaussian(x,-3,1.2,1.0)
    out=run_multistate_exact(psi,x,dx,V,mass=20,dt=0.002,steps=100,store_every=10)
    assert np.max(np.abs(out["norm"]-1)) < 1e-11


def test_gauge_alignment_and_nac_antisymmetry():
    x,dx=uniform_grid(-6,6,401)
    E,U,g,d=adiabatic_grid(x)
    overlaps=np.einsum("xki,xki->xi", U[:-1], U[1:])
    assert np.min(overlaps) > 0
    assert np.max(np.abs(d[:,0,1]+d[:,1,0])) < 1e-12


def test_covariant_adiabatic_operator_matches_transformed_diabatic_interior():
    x,dx=uniform_grid(-8,8,2401)
    V=avoided_crossing_diabatic(x)
    E,U,g,d=adiabatic_grid(x)

    chi=np.zeros((len(x),2),complex)
    envelope=np.exp(-0.25*x*x)
    chi[:,0]=envelope*np.exp(0.35j*x)
    chi[:,1]=0.2*envelope*np.exp(-0.2j*x)

    Hd_ad=adiabatic_hamiltonian_action(chi,x,dx,5.0,E,d)

    psi=transform_ad_to_dia(chi,U)
    Hd_d=diabatic_hamiltonian_action(psi,x,dx,5.0,V)
    transformed=transform_dia_to_ad(Hd_d,U)

    sl=slice(20,-20)
    err=np.sqrt(np.mean(np.abs(Hd_ad[sl]-transformed[sl])**2))
    ref=np.sqrt(np.mean(np.abs(transformed[sl])**2))
    assert err/ref < 2e-3


def test_energy_conserving_spawn_when_allowed():
    parent=AdiabaticTBF(state=1,q=0.0,p=1.0,alpha=1.0)
    child=spawn_child_energy_conserving(parent,mass=20.0)
    assert child is not None

    from gaussian_dynamics.adiabatic import adiabatic_point
    E,_,_,_=adiabatic_point(0.0)
    before=parent.p**2/(40.0)+E[parent.state]
    after=child.p**2/(40.0)+E[child.state]
    assert abs(before-after) < 1e-12


def test_spawn_basis_growth_preserves_old_coefficients():
    old=[AdiabaticTBF(state=1,q=0.0,p=1.0,alpha=1.0)]
    new=maybe_spawn(old,threshold=1e-6,mass=20.0)
    assert len(new) >= len(old)
    C=expand_coefficients_for_new_basis(np.array([1+0j]),old,new)
    assert C[0] == 1+0j
    if len(new)>1:
        assert np.allclose(C[1:],0)


def test_adiabatic_gaussian_H_is_hermitian_to_quadrature_accuracy():
    x,dx=uniform_grid(-8,8,1401)
    basis=[
        AdiabaticTBF(0,-0.5,0.4,1.0),
        AdiabaticTBF(1,0.2,0.1,1.0),
    ]
    S,H=adiabatic_gaussian_basis_matrices(x,dx,basis,mass=10.0)
    assert np.allclose(S,S.conj().T,atol=1e-10)
    assert np.allclose(H,H.conj().T,atol=2e-4)
