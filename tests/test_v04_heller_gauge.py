import numpy as np

from gaussian_dynamics.heller_nd import run_thawed_gaussian_nd
from gaussian_dynamics.gauge_transport import align_subspace, projector_distance


def test_2d_harmonic_coherent_width_is_constant_and_center_exact():
    mass=2.0
    omega=0.7
    sigma=np.sqrt(1.0/(2.0*mass*omega))
    q0=np.array([-0.8,0.5])
    p0=np.array([0.4,-0.3])

    V=lambda q: 0.5*mass*omega**2*np.dot(q,q)
    G=lambda q: mass*omega**2*np.asarray(q)
    H=lambda q: mass*omega**2*np.eye(2)

    dt=0.001
    steps=700
    t=dt*steps

    out=run_thawed_gaussian_nd(
        q0,p0,sigma,mass,V,G,H,dt=dt,steps=steps,store_every=steps
    )

    q_exact=q0*np.cos(omega*t)+p0/(mass*omega)*np.sin(omega*t)
    A_expected=1j*mass*omega*np.eye(2)

    assert np.linalg.norm(out["q"][-1]-q_exact) < 2e-11
    assert np.linalg.norm(out["A"][-1]-A_expected) < 2e-11
    assert np.linalg.norm(out["A"][-1]-out["A"][-1].T) < 1e-14


def test_subspace_procrustes_removes_arbitrary_unitary_rotation():
    rng=np.random.default_rng(7)
    X=rng.normal(size=(5,2))+1j*rng.normal(size=(5,2))
    Q,_=np.linalg.qr(X)

    Z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    W,_=np.linalg.qr(Z)

    current=Q@W
    aligned,_=align_subspace(Q,current)

    assert projector_distance(Q,current) < 1e-12
    assert np.linalg.norm(aligned-Q) < 1e-12
