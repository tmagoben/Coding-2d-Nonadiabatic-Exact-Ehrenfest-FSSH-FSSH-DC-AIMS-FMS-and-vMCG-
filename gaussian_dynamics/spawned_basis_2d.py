from dataclasses import dataclass
import numpy as np

from .ci2d import (
    LVC2DParameters,
    adiabatic_energies_2d,
    adiabatic_gradients_2d,
    vector_nac_2d,
)
from .gaussian_nd import (
    gaussian_nd,
    gaussian_nd_gradient,
    gaussian_nd_laplacian,
    gaussian_nd_time_derivative,
    analytic_overlap_equal_width,
)


@dataclass
class AdiabaticTBF2D:
    state: int
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray

    def copy(self):
        return AdiabaticTBF2D(
            int(self.state),
            np.asarray(self.q,float).copy(),
            np.asarray(self.p,float).copy(),
            np.asarray(self.A,float).copy(),
        )


def midpoint_grid_2d(xmin=-5.0, xmax=5.0, nx=48, ymin=-5.0, ymax=5.0, ny=48):
    """Cell-centered grid. For even nx,ny the exact CI at (0,0) is excluded."""
    dx=(xmax-xmin)/nx
    dy=(ymax-ymin)/ny
    x=xmin+(np.arange(nx)+0.5)*dx
    y=ymin+(np.arange(ny)+0.5)*dy
    X,Y=np.meshgrid(x,y,indexing="ij")
    points=np.stack([X,Y],axis=-1)
    return x,y,X,Y,points,dx,dy


def adiabatic_fields_2d(X, Y, p=LVC2DParameters()):
    """Return E[...,state], d[...,i,j,coord], tau[...,i,j].

    tau = div(d) + sum_alpha d_alpha d_alpha, evaluated numerically.
    The grid must avoid the exact CI.
    """
    E = adiabatic_energies_2d(X,Y,p)

    denom=(p.kappa*X)**2+(p.lam*Y)**2
    if np.any(denom <= 0.0):
        raise ValueError("Grid contains the exact CI; use a cell-centered grid.")

    a_x = -0.5*p.kappa*p.lam*Y/denom
    a_y =  0.5*p.kappa*p.lam*X/denom

    d=np.zeros(X.shape+(2,2,2),float)
    d[...,0,1,0]=a_x
    d[...,0,1,1]=a_y
    d[...,1,0,0]=-a_x
    d[...,1,0,1]=-a_y

    # infer spacings from mesh
    dx=float(X[1,0]-X[0,0])
    dy=float(Y[0,1]-Y[0,0])

    ddx=np.gradient(d[...,0],dx,axis=0,edge_order=2)
    ddy=np.gradient(d[...,1],dy,axis=1,edge_order=2)
    div=ddx+ddy

    d2=(
        np.einsum("...ik,...kj->...ij",d[...,0],d[...,0])
        +np.einsum("...ik,...kj->...ij",d[...,1],d[...,1])
    )
    tau=div+d2
    return E,d,tau


def tbf_guidance(tbf, mass, p=LVC2DParameters()):
    grad=adiabatic_gradients_2d(tbf.q,p)
    qdot=np.asarray(tbf.p,float)/mass
    pdot=-grad[tbf.state]
    return qdot,pdot


def basis_scalar_functions(points,basis):
    return [
        gaussian_nd(points,b.q,b.p,b.A)
        for b in basis
    ]


def basis_matrices_2d(points, dx, dy, basis, mass=20.0, p=LVC2DParameters()):
    """Build S,H,T for adiabatic TBFs with the covariant adiabatic Hamiltonian."""
    X=points[...,0]
    Y=points[...,1]
    E,d,tau_el=adiabatic_fields_2d(X,Y,p)

    n=len(basis)
    S=np.zeros((n,n),complex)
    H=np.zeros((n,n),complex)
    T=np.zeros((n,n),complex)

    gs=[]; grads=[]; laps=[]; gdots=[]
    for b in basis:
        g=gaussian_nd(points,b.q,b.p,b.A)
        gg=gaussian_nd_gradient(points,b.q,b.p,b.A)
        lap=gaussian_nd_laplacian(points,b.q,b.p,b.A)
        qdot,pdot=tbf_guidance(b,mass,p)
        gdot=gaussian_nd_time_derivative(points,b.q,b.p,b.A,qdot,pdot)

        gs.append(g); grads.append(gg); laps.append(lap); gdots.append(gdot)

    H_on_basis=[]
    for j,b in enumerate(basis):
        out=np.zeros(points.shape[:-1]+(2,),complex)
        for a in (0,1):
            first = 2.0*(
                d[...,a,b.state,0]*grads[j][...,0]
                +d[...,a,b.state,1]*grads[j][...,1]
            )
            second=tau_el[...,a,b.state]*gs[j]
            base=laps[j] if a==b.state else 0.0
            out[...,a]=-(base+first+second)/(2.0*mass)
            if a==b.state:
                out[...,a]+=E[...,a]*gs[j]
        H_on_basis.append(out)

    weight=dx*dy
    for i,bi in enumerate(basis):
        for j,bj in enumerate(basis):
            if bi.state==bj.state:
                S[i,j]=np.vdot(gs[i],gs[j])*weight
                T[i,j]=np.vdot(gs[i],gdots[j])*weight

            # electronic inner product selects component bi.state
            H[i,j]=np.vdot(gs[i],H_on_basis[j][...,bi.state])*weight

    return S,H,T


def coefficient_rhs(C,S,H,T):
    rhs=-1j*(H@C)-T@C
    return np.linalg.solve(S,rhs)


def nac_coupling_indicator(tbf,mass,p=LVC2DParameters()):
    target=1-tbf.state
    d=vector_nac_2d(tbf.q,p)[tbf.state,target]
    return float(abs(np.dot(np.asarray(tbf.p)/mass,d)))


def energy_conserving_child_nac(tbf,mass,p=LVC2DParameters()):
    target=1-tbf.state
    d=vector_nac_2d(tbf.q,p)[tbf.state,target]
    dn=np.linalg.norm(d)
    if dn < 1e-14:
        return None

    n=d/dn
    E=adiabatic_energies_2d(tbf.q[0],tbf.q[1],p)
    delta=E[target]-E[tbf.state]
    pn=float(np.dot(tbf.p,n))
    disc=pn*pn-2.0*mass*delta
    if disc < 0:
        return None

    root=np.sqrt(disc)
    candidates=[-pn+root,-pn-root]
    eta=min(candidates,key=abs)
    pchild=np.asarray(tbf.p,float)+eta*n

    return AdiabaticTBF2D(
        target,
        np.asarray(tbf.q,float).copy(),
        pchild,
        np.asarray(tbf.A,float).copy(),
    )


def phase_space_overlap_magnitude(a,b):
    if a.state != b.state:
        return 0.0
    if not np.allclose(a.A,b.A,atol=1e-12):
        return 0.0
    return float(abs(analytic_overlap_equal_width(a.q,a.p,b.q,b.p,a.A)))


def maybe_spawn_once(
    basis,
    threshold,
    mass,
    overlap_block=0.85,
    p=LVC2DParameters(),
):
    """At most one deterministic spawn per call."""
    for idx,parent in enumerate(basis):
        if nac_coupling_indicator(parent,mass,p) <= threshold:
            continue

        child=energy_conserving_child_nac(parent,mass,p)
        if child is None:
            continue

        redundant=any(
            b.state==child.state
            and phase_space_overlap_magnitude(b,child)>=overlap_block
            for b in basis
        )
        if not redundant:
            return idx, child
    return None, None


def _temporary_basis(states, As, q, mom):
    return [
        AdiabaticTBF2D(int(states[i]),q[i].copy(),mom[i].copy(),As[i].copy())
        for i in range(len(states))
    ]


def run_coupled_spawned_basis_2d(
    points,
    dx,
    dy,
    initial_basis,
    C0,
    mass=20.0,
    dt=0.0005,
    steps=100,
    spawn_threshold=0.02,
    overlap_block=0.85,
    max_basis=6,
    p=LVC2DParameters(),
    store_every=5,
):
    """Coupled moving-basis propagation with dynamic zero-amplitude spawning.

    This is a transparent FMS/AIMS-style prototype on the analytic 2D CI model.
    """
    basis=[b.copy() for b in initial_basis]
    C=np.asarray(C0,complex).copy()

    S,_,_=basis_matrices_2d(points,dx,dy,basis,mass,p)
    norm0=np.real(np.vdot(C,S@C))
    C/=np.sqrt(norm0)

    times=[]; norms=[]; sizes=[]; pops=[]; conds=[]; events=[]

    def diagnostics(step):
        S,H,T=basis_matrices_2d(points,dx,dy,basis,mass,p)
        norm=np.real(np.vdot(C,S@C))
        state_pop=np.zeros(2)
        for state in (0,1):
            idx=[i for i,b in enumerate(basis) if b.state==state]
            if idx:
                block=S[np.ix_(idx,idx)]
                cc=C[idx]
                state_pop[state]=np.real(np.vdot(cc,block@cc))
        if norm>0:
            state_pop/=norm

        times.append(step*dt)
        norms.append(norm)
        sizes.append(len(basis))
        pops.append(state_pop)
        conds.append(np.linalg.cond(S))

    diagnostics(0)

    for step in range(1,steps+1):
        n=len(basis)
        states=np.array([b.state for b in basis],int)
        As=[b.A.copy() for b in basis]
        q=np.array([b.q for b in basis],float)
        mom=np.array([b.p for b in basis],float)

        def rhs(Cx,qx,px):
            temp=_temporary_basis(states,As,qx,px)
            S,H,T=basis_matrices_2d(points,dx,dy,temp,mass,p)
            Cdot=coefficient_rhs(Cx,S,H,T)
            qdot=np.zeros_like(qx)
            pdot=np.zeros_like(px)
            for i,b in enumerate(temp):
                qdot[i],pdot[i]=tbf_guidance(b,mass,p)
            return Cdot,qdot,pdot

        k1=rhs(C,q,mom)
        k2=rhs(C+0.5*dt*k1[0],q+0.5*dt*k1[1],mom+0.5*dt*k1[2])
        k3=rhs(C+0.5*dt*k2[0],q+0.5*dt*k2[1],mom+0.5*dt*k2[2])
        k4=rhs(C+dt*k3[0],q+dt*k3[1],mom+dt*k3[2])

        C=C+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6.0
        q=q+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6.0
        mom=mom+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6.0

        for i,b in enumerate(basis):
            b.q=q[i].copy()
            b.p=mom[i].copy()

        if len(basis) < max_basis:
            parent_idx,child=maybe_spawn_once(
                basis,spawn_threshold,mass,overlap_block,p
            )
            if child is not None:
                basis.append(child)
                C=np.concatenate([C,[0.0+0.0j]])
                events.append({
                    "step":step,
                    "time":step*dt,
                    "parent_index":int(parent_idx),
                    "new_index":len(basis)-1,
                    "target_state":int(child.state),
                })

        if step % store_every == 0:
            diagnostics(step)

    return {
        "time":np.asarray(times),
        "norm":np.asarray(norms),
        "basis_size":np.asarray(sizes),
        "state_populations":np.asarray(pops),
        "condition_number":np.asarray(conds),
        "events":events,
        "final_basis":basis,
        "final_coefficients":C,
    }
