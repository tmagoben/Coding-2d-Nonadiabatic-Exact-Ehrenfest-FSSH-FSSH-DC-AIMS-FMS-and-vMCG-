import numpy as np
from .gauge_graph import nearest_unitary

def hermitian_exponential(H, dt):
    H=np.asarray(H,complex)
    if H.ndim!=2 or H.shape[0]!=H.shape[1]: raise ValueError('H must be square')
    if not np.allclose(H,H.conj().T,atol=1e-10): raise ValueError('H must be Hermitian')
    e,U=np.linalg.eigh(H)
    return (U*np.exp(-1j*dt*e))@U.conj().T

def explicit_nac_step(coefficients, energies, nac, velocity, dt):
    c=np.asarray(coefficients,complex); E=np.asarray(energies,float)
    d=np.asarray(nac,float); v=np.asarray(velocity,float); ns=len(E)
    if c.shape!=(ns,): raise ValueError('coefficient vector has incompatible dimension')
    if d.ndim!=3 or d.shape[:2]!=(ns,ns): raise ValueError('nac must have shape (nstate,nstate,nq)')
    if d.shape[2]!=len(v): raise ValueError('velocity and NAC dimensions differ')
    directional=np.einsum('ija,a->ij',d,v)
    Heff=np.diag(E).astype(complex)-1j*directional
    return hermitian_exponential(Heff,dt)@c

def overlap_strang_step(coefficients,H_old,H_new,overlap_old_new,dt):
    c=np.asarray(coefficients,complex); H0=np.asarray(H_old,complex); H1=np.asarray(H_new,complex); O=np.asarray(overlap_old_new,complex)
    if H0.shape!=H1.shape or H0.shape[0]!=H0.shape[1]: raise ValueError('H_old and H_new must be equal-size square matrices')
    if c.shape!=(H0.shape[0],): raise ValueError('coefficient vector has incompatible dimension')
    if O.shape!=H0.shape: raise ValueError('overlap has incompatible dimension')
    U=nearest_unitary(O)
    c=hermitian_exponential(H0,0.5*dt)@c
    c=U.conj().T@c
    c=hermitian_exponential(H1,0.5*dt)@c
    return c

def overlap_step_operator(H_old,H_new,overlap_old_new,dt):
    H0=np.asarray(H_old,complex); H1=np.asarray(H_new,complex); U=nearest_unitary(np.asarray(overlap_old_new,complex))
    return hermitian_exponential(H1,0.5*dt)@U.conj().T@hermitian_exponential(H0,0.5*dt)

def electronic_fidelity(c1,c2):
    a=np.asarray(c1,complex); b=np.asarray(c2,complex); na=np.vdot(a,a).real; nb=np.vdot(b,b).real
    if na<=0 or nb<=0: raise ValueError('zero electronic vector')
    return float(abs(np.vdot(a,b))**2/(na*nb))
