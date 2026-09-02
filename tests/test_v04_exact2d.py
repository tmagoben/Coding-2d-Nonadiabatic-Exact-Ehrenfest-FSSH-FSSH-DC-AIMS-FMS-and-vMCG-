import numpy as np

from gaussian_dynamics.ci2d import diabatic_potential_2d
from gaussian_dynamics.gaussian_nd import gaussian_nd
from gaussian_dynamics.exact2d import run_exact_2d


def test_exact_2d_two_state_norm_conservation():
    n=40
    L=5.0
    dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    y=x.copy()
    X,Y=np.meshgrid(x,y,indexing="ij")
    P=np.stack([X,Y],axis=-1)

    g=gaussian_nd(
        P,
        q=np.array([-2.0,0.6]),
        p=np.array([1.0,-0.1]),
        A=np.eye(2),
    )

    psi=np.zeros((n,n,2),complex)
    psi[...,1]=g
    V=diabatic_potential_2d(X,Y)

    out=run_exact_2d(
        psi,dx,dx,V,mass=20.0,dt=0.004,steps=60,store_every=10
    )

    assert np.max(np.abs(out["norm"]-1.0)) < 2e-12
    assert np.max(np.abs(np.sum(out["populations"],axis=1)-1.0)) < 2e-12


def test_exact_2d_strang_self_convergence():
    n=24
    L=4.0
    dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    X,Y=np.meshgrid(x,x,indexing="ij")
    P=np.stack([X,Y],axis=-1)

    g=gaussian_nd(P,np.array([-1.5,0.5]),np.array([0.8,-0.1]),np.eye(2))
    psi=np.zeros((n,n,2),complex)
    psi[...,1]=g
    V=diabatic_potential_2d(X,Y)

    finals=[]
    for dt,steps in [(0.008,30),(0.004,60),(0.002,120)]:
        out=run_exact_2d(
            psi,dx,dx,V,mass=20.0,dt=dt,steps=steps,store_every=steps
        )
        finals.append(out["psi"][-1])

    def phase_aligned_error(a,ref):
        inner=np.sum(np.conj(ref)*a)*dx*dx
        phase=inner/abs(inner)
        return np.sqrt(np.sum(np.abs(a-phase*ref)**2)*dx*dx)

    coarse=phase_aligned_error(finals[0],finals[2])
    medium=phase_aligned_error(finals[1],finals[2])

    # With the finest run used as a numerical reference, second-order Strang
    # convergence should make the coarser error substantially larger.
    assert coarse/medium > 3.0
