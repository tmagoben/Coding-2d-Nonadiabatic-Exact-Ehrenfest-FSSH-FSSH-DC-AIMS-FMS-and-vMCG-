import numpy as np
from gaussian_dynamics.ci2d import adiabatic_energies_2d,analytic_adiabatic_vectors,vector_nac_2d
from gaussian_dynamics.temporal_electronic import explicit_nac_step,overlap_strang_step,electronic_fidelity

def random_unitary(rng,n):
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));Q,_=np.linalg.qr(X);return Q

def test_overlap_strang_is_locally_gauge_covariant():
    rng=np.random.default_rng(5);H0=np.array([[.1,.02],[.02,.3]],complex);H1=np.array([[.12,-.01],[-.01,.28]],complex);O=random_unitary(rng,2);c=np.array([.8+.1j,.2-.3j]);c/=np.linalg.norm(c)
    out=overlap_strang_step(c,H0,H1,O,.03);G0=random_unitary(rng,2);G1=random_unitary(rng,2)
    outp=overlap_strang_step(G0.conj().T@c,G0.conj().T@H0@G0,G1.conj().T@H1@G1,G0.conj().T@O@G1,.03)
    assert np.allclose(outp,G1.conj().T@out,atol=2e-12)

def prop(nsteps):
    q0=np.array([-.8,.65]);q1=np.array([.8,.65]);T=.7;dt=T/nsteps;v=(q1-q0)/T;cn=np.array([1+0j,0j]);co=cn.copy()
    for n in range(nsteps):
        qa=q0+n/nsteps*(q1-q0);qb=q0+(n+1)/nsteps*(q1-q0);qm=.5*(qa+qb)
        cn=explicit_nac_step(cn,adiabatic_energies_2d(*qm),vector_nac_2d(qm),v,dt)
        Ua=analytic_adiabatic_vectors(qa);Ub=analytic_adiabatic_vectors(qb)
        co=overlap_strang_step(co,np.diag(adiabatic_energies_2d(*qa)),np.diag(adiabatic_energies_2d(*qb)),Ua.conj().T@Ub,dt)
    return cn,co

def test_overlap_and_explicit_nac_converge_to_same_result():
    c1,o1=prop(100);c2,o2=prop(400)
    assert electronic_fidelity(c1,c2)>.999999
    assert electronic_fidelity(o1,o2)>.999999
    assert electronic_fidelity(c2,o2)>.99999
    assert abs(np.linalg.norm(c2)-1)<1e-12 and abs(np.linalg.norm(o2)-1)<1e-12
