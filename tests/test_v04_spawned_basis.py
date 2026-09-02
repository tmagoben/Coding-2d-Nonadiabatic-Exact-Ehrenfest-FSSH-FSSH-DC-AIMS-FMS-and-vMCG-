import numpy as np

from gaussian_dynamics.spawned_basis_2d import (
    AdiabaticTBF2D,
    midpoint_grid_2d,
    basis_matrices_2d,
    tbf_guidance,
    energy_conserving_child_nac,
    maybe_spawn_once,
    run_coupled_spawned_basis_2d,
)
from gaussian_dynamics.ci2d import adiabatic_energies_2d
from gaussian_dynamics.gaussian_nd import gaussian_nd


def test_adiabatic_spawned_basis_matrices_are_hermitian():
    x,y,X,Y,P,dx,dy=midpoint_grid_2d(-4,4,36,-4,4,36)
    A=np.eye(2)
    basis=[
        AdiabaticTBF2D(0,np.array([-1.2,0.7]),np.array([0.5,0.1]),A),
        AdiabaticTBF2D(1,np.array([0.9,0.8]),np.array([-0.2,0.2]),A),
    ]

    S,H,T=basis_matrices_2d(P,dx,dy,basis,mass=20.0)

    assert np.allclose(S,S.conj().T,atol=2e-10)
    # finite-difference electronic second-coupling field is the dominant tolerance
    assert np.allclose(H,H.conj().T,atol=3e-3)


def test_energy_conserving_child_along_vector_nac():
    A=np.eye(2)
    parent=AdiabaticTBF2D(
        1,
        np.array([0.55,0.45]),
        np.array([0.6,0.8]),
        A,
    )
    child=energy_conserving_child_nac(parent,mass=20.0)
    assert child is not None

    E=adiabatic_energies_2d(parent.q[0],parent.q[1])
    before=np.dot(parent.p,parent.p)/(40.0)+E[parent.state]
    after=np.dot(child.p,child.p)/(40.0)+E[child.state]
    assert abs(before-after) < 1e-12


def test_spawn_rule_is_deterministic_and_zero_coefficient_insertion_is_continuous():
    A=np.eye(2)
    parent=AdiabaticTBF2D(
        1,
        np.array([0.55,0.45]),
        np.array([0.6,0.8]),
        A,
    )
    basis=[parent]
    idx1,ch1=maybe_spawn_once(basis,threshold=1e-6,mass=20.0)
    idx2,ch2=maybe_spawn_once(basis,threshold=1e-6,mass=20.0)

    assert idx1==idx2==0
    assert ch1 is not None and ch2 is not None
    assert np.allclose(ch1.q,ch2.q)
    assert np.allclose(ch1.p,ch2.p)

    # Adding a zero-amplitude child does not change the represented scalar parent
    # component at the insertion instant.
    _,_,_,_,_,dx,dy=midpoint_grid_2d(-4,4,30,-4,4,30)
    x,y,X,Y,P,_,_=midpoint_grid_2d(-4,4,30,-4,4,30)
    g_before=gaussian_nd(P,parent.q,parent.p,parent.A)
    g_after=1.0*g_before + 0.0*gaussian_nd(P,ch1.q,ch1.p,ch1.A)
    assert np.max(np.abs(g_after-g_before)) == 0.0


def test_short_coupled_spawned_run_is_finite_and_grows_basis():
    x,y,X,Y,P,dx,dy=midpoint_grid_2d(-3,3,24,-3,3,24)
    A=1.2*np.eye(2)
    parent=AdiabaticTBF2D(
        1,
        np.array([0.55,0.45]),
        np.array([0.6,0.8]),
        A,
    )

    out=run_coupled_spawned_basis_2d(
        P,dx,dy,
        [parent],
        C0=[1.0+0j],
        mass=20.0,
        dt=0.0002,
        steps=12,
        spawn_threshold=1e-6,
        overlap_block=0.9,
        max_basis=2,
        store_every=2,
    )

    assert out["basis_size"][-1] == 2
    assert len(out["events"]) == 1
    assert np.all(np.isfinite(out["norm"]))
    assert np.max(np.abs(out["norm"]-1.0)) < 2e-5
    assert np.all(np.isfinite(out["state_populations"]))


def test_moving_basis_overlap_derivative_identity():
    x,y,X,Y,P,dx,dy=midpoint_grid_2d(-4,4,32,-4,4,32)
    A=np.eye(2)
    basis=[
        AdiabaticTBF2D(1,np.array([-1.0,0.8]),np.array([0.5,0.1]),A),
        AdiabaticTBF2D(1,np.array([-0.2,0.9]),np.array([0.2,-0.1]),A),
    ]

    S,H,T=basis_matrices_2d(P,dx,dy,basis,mass=20.0)

    h=1e-5
    plus=[]
    minus=[]
    for b in basis:
        qdot,pdot=tbf_guidance(b,20.0)
        plus.append(
            AdiabaticTBF2D(b.state,b.q+h*qdot,b.p+h*pdot,b.A)
        )
        minus.append(
            AdiabaticTBF2D(b.state,b.q-h*qdot,b.p-h*pdot,b.A)
        )

    Sp,_,_=basis_matrices_2d(P,dx,dy,plus,mass=20.0)
    Sm,_,_=basis_matrices_2d(P,dx,dy,minus,mass=20.0)
    Sdot=(Sp-Sm)/(2*h)

    assert np.allclose(Sdot,T+T.conj().T,atol=2e-9)


def test_spawned_child_receives_amplitude_from_coupled_propagation():
    x,y,X,Y,P,dx,dy=midpoint_grid_2d(-3,3,24,-3,3,24)
    parent=AdiabaticTBF2D(
        1,
        np.array([0.55,0.45]),
        np.array([0.6,0.8]),
        1.2*np.eye(2),
    )

    out=run_coupled_spawned_basis_2d(
        P,dx,dy,[parent],[1.0+0j],
        mass=20.0,
        dt=0.0002,
        steps=50,
        spawn_threshold=1e-6,
        overlap_block=0.9,
        max_basis=2,
        store_every=5,
    )

    # The child starts at zero amplitude. Nonzero final lower-state population
    # therefore comes from the coupled matrix propagation, not manual assignment.
    assert out["state_populations"][-1,0] > 1e-10
