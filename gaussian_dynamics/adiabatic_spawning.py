from dataclasses import dataclass
import numpy as np

from .gaussian import frozen_gaussian, kinetic_on_gaussian
from .grids import inner_product
from .adiabatic import adiabatic_point, adiabatic_grid, adiabatic_hamiltonian_action


@dataclass(frozen=True)
class AdiabaticTBF:
    state: int
    q: float
    p: float
    alpha: float = 1.0


def coupling_indicator(tbf, mass=1.0):
    E,U,g,d = adiabatic_point(tbf.q)
    target = 1-tbf.state
    return abs((tbf.p/mass)*d[tbf.state,target])


def spawn_child_energy_conserving(parent, mass=1.0):
    target = 1-parent.state
    E,_,_,_ = adiabatic_point(parent.q)
    rad = parent.p**2 + 2.0*mass*(E[parent.state]-E[target])
    if rad < 0:
        return None
    sign = 1.0 if parent.p >= 0 else -1.0
    return AdiabaticTBF(target, parent.q, sign*np.sqrt(rad), parent.alpha)


def tbf_overlap_same_state(a, b):
    if a.state != b.state:
        return 0.0
    dq=a.q-b.q; dp=a.p-b.p; alpha=a.alpha
    if abs(a.alpha-b.alpha) > 1e-12:
        return 0.0
    return float(np.exp(-0.25*alpha*dq*dq - 0.25*dp*dp/alpha))


def maybe_spawn(basis, threshold=0.005, overlap_block=0.85, mass=1.0):
    """Deterministic one-pass spawning rule. Returns a new list."""
    new_basis = list(basis)
    for parent in list(basis):
        if coupling_indicator(parent,mass) <= threshold:
            continue
        child = spawn_child_energy_conserving(parent,mass)
        if child is None:
            continue
        blocked = any(
            b.state == child.state and tbf_overlap_same_state(b,child) >= overlap_block
            for b in new_basis
        )
        if not blocked:
            new_basis.append(child)
    return new_basis


def adiabatic_gaussian_basis_matrices(x, dx, basis, mass=1.0):
    """S and H for basis g_k(x)|phi_state(x)> using exact adiabatic kinetic couplings."""
    E,U,grad,d = adiabatic_grid(x)
    n=len(basis)
    S=np.zeros((n,n),complex)
    H=np.zeros((n,n),complex)

    basis_vectors=[]
    for b in basis:
        vec=np.zeros((len(x),2),complex)
        vec[:,b.state]=frozen_gaussian(x,b.q,b.p,b.alpha)
        basis_vectors.append(vec)

    Hvec=[adiabatic_hamiltonian_action(v,x,dx,mass,E,d) for v in basis_vectors]

    for i,bi in enumerate(basis):
        for j,bj in enumerate(basis):
            S[i,j]=np.sum(np.conj(basis_vectors[i])*basis_vectors[j])*dx
            H[i,j]=np.sum(np.conj(basis_vectors[i])*Hvec[j])*dx
    return S,H


def expand_coefficients_for_new_basis(C_old, old_basis, new_basis):
    """Newly spawned functions receive zero amplitude, preserving the old wavefunction."""
    C=np.zeros(len(new_basis),complex)
    C[:len(old_basis)] = C_old
    return C
