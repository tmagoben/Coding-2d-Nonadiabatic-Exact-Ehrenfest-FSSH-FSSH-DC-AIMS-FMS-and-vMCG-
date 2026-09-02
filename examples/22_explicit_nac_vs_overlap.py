import numpy as np
from gaussian_dynamics.ci2d import adiabatic_energies_2d,analytic_adiabatic_vectors,vector_nac_2d
from gaussian_dynamics.temporal_electronic import explicit_nac_step,overlap_strang_step,electronic_fidelity

q0=np.array([-0.8,0.65]); q1=np.array([0.8,0.65]); total_time=0.7; nsteps=400; dt=total_time/nsteps; velocity=(q1-q0)/total_time
c_nac=np.array([1.0+0j,0.0+0j]); c_overlap=c_nac.copy()
for n in range(nsteps):
    qa=q0+n/nsteps*(q1-q0); qb=q0+(n+1)/nsteps*(q1-q0); qm=0.5*(qa+qb)
    c_nac=explicit_nac_step(c_nac,adiabatic_energies_2d(*qm),vector_nac_2d(qm),velocity,dt)
    Ua=analytic_adiabatic_vectors(qa); Ub=analytic_adiabatic_vectors(qb)
    c_overlap=overlap_strang_step(c_overlap,np.diag(adiabatic_energies_2d(*qa)),np.diag(adiabatic_energies_2d(*qb)),Ua.conj().T@Ub,dt)
print('Explicit-NAC coefficients:',c_nac)
print('Overlap/local-diabatic coefficients:',c_overlap)
print('Fidelity:',electronic_fidelity(c_nac,c_overlap))
print('Norms:',np.linalg.norm(c_nac),np.linalg.norm(c_overlap))
