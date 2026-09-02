import numpy as np
from .local_gaussian_nd import basis_time_matrix_element_equal_width

def nuclear_seed_basis_time_matrix(basis,registry,reference_selector,qdots,pdots):
    n=len(basis); qdots=np.asarray(qdots,float); pdots=np.asarray(pdots,float)
    if qdots.shape!=(n,len(basis[0].q)) or pdots.shape!=qdots.shape: raise ValueError('kinematic arrays have incompatible shape')
    T=np.zeros((n,n),complex)
    for i in range(n):
        for j in range(n):
            if not np.allclose(basis[i].A,basis[j].A,atol=1e-12): raise ValueError('equal widths required')
            ref=reference_selector(i,j)
            fac=registry.pair_factors(basis[i].node,basis[i].electronic_coefficients,basis[j].node,basis[j].electronic_coefficients,ref)
            Tn=basis_time_matrix_element_equal_width(basis[i].q,basis[i].p,basis[j].q,basis[j].p,basis[j].A,qdots[j],pdots[j])
            T[i,j]=Tn*fac['overlap']
    return T

def metric_compatible_basis_connection(S_old,S_new,dt,seed=None):
    S0=np.asarray(S_old,complex); S1=np.asarray(S_new,complex)
    if S0.shape!=S1.shape or S0.shape[0]!=S0.shape[1]: raise ValueError('S matrices incompatible')
    if dt<=0: raise ValueError('dt must be positive')
    Sdot=(S1-S0)/dt
    if seed is None: return 0.5*Sdot
    seed=np.asarray(seed,complex)
    if seed.shape!=S0.shape: raise ValueError('seed incompatible')
    return seed+0.5*(Sdot-seed-seed.conj().T)

def basis_connection_residual(S_old,S_new,T,dt):
    Sdot=(np.asarray(S_new)-np.asarray(S_old))/dt; T=np.asarray(T)
    return float(np.linalg.norm(Sdot-T-T.conj().T,ord='fro'))

def linear_rk4_step(C,generator,dt):
    C=np.asarray(C,complex); A=np.asarray(generator,complex)
    k1=A@C; k2=A@(C+0.5*dt*k1); k3=A@(C+0.5*dt*k2); k4=A@(C+dt*k3)
    return C+dt*(k1+2*k2+2*k3+k4)/6

def moving_basis_coefficient_step(C,S_old,H_old,S_new,H_new,T_mid,dt):
    Sm=0.5*(np.asarray(S_old)+np.asarray(S_new)); Hm=0.5*(np.asarray(H_old)+np.asarray(H_new)); T=np.asarray(T_mid)
    A=np.linalg.solve(Sm,-1j*Hm-T); return linear_rk4_step(C,A,dt)

def generalized_norm(C,S):
    C=np.asarray(C,complex); S=np.asarray(S,complex); return float(np.real(np.vdot(C,S@C)))
