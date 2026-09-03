"""Independent moving-frame validation oracles for v0.28.0."""
from dataclasses import dataclass
import itertools
import numpy as np
from scipy.linalg import expm
from .moving_frame_v280 import moving_frame_hamiltonian_v280


def _scaled_norm_v280(left, right):
    left=np.asarray(left); right=np.asarray(right)
    if left.shape != right.shape: return float('inf')
    return float(np.linalg.norm(left-right)/max(float(np.linalg.norm(left)),float(np.linalg.norm(right)),1.0))


def periodic_second_derivative_v280(count, spacing):
    count=int(count); spacing=float(spacing)
    if count < 3 or not np.isfinite(spacing) or spacing <= 0: raise ValueError('periodic lattice requires at least three points and positive spacing.')
    m=np.zeros((count,count),float)
    for i in range(count):
        m[i,i]=-2/spacing**2; m[i,(i-1)%count]=1/spacing**2; m[i,(i+1)%count]=1/spacing**2
    return m


def scalar_lattice_kinetic_v280(axes, mass_matrix_au):
    axes=tuple(np.asarray(a,float) for a in axes); ndim=len(axes); mass=np.asarray(mass_matrix_au,float)
    if mass.shape != (ndim,ndim) or not np.allclose(mass,np.diag(np.diag(mass)),atol=1e-13,rtol=0):
        raise ValueError('independent v0.28 lattice oracle currently requires a diagonal mass matrix.')
    shapes=tuple(len(a) for a in axes)
    if any(n<3 for n in shapes): raise ValueError('every periodic lattice axis requires at least three points.')
    kinetic=np.zeros((int(np.prod(shapes)),int(np.prod(shapes))))
    for k,(axis,mval) in enumerate(zip(axes,np.diag(mass))):
        diffs=np.diff(axis)
        if len(diffs)==0 or not np.allclose(diffs,diffs[0],atol=1e-13,rtol=1e-12): raise ValueError('lattice axes must be uniformly spaced.')
        d2=periodic_second_derivative_v280(len(axis),diffs[0]); factors=[d2 if j==k else np.eye(n) for j,n in enumerate(shapes)]
        block=factors[0]
        for factor in factors[1:]: block=np.kron(block,factor)
        kinetic += -0.5*block/float(mval)
    return kinetic


def lattice_points_v280(axes):
    axes=tuple(np.asarray(a,float) for a in axes); mesh=np.meshgrid(*axes,indexing='ij')
    return np.stack([x.reshape(-1) for x in mesh],axis=-1)


@dataclass(frozen=True)
class LatticeGaugeOracleV280:
    points: np.ndarray
    fixed_hamiltonian: np.ndarray
    moving_hamiltonian: np.ndarray
    similarity_hamiltonian: np.ndarray
    transformation: np.ndarray
    maximum_link_unitarity_residual: float
    hamiltonian_hermiticity_residual: float
    similarity_residual: float
    def validate(self,tolerance=3e-12):
        for m in (self.fixed_hamiltonian,self.moving_hamiltonian,self.similarity_hamiltonian,self.transformation):
            if m.ndim!=2 or m.shape[0]!=m.shape[1] or not np.all(np.isfinite(m)): raise ValueError('lattice oracle contains an invalid matrix.')
        if self.fixed_hamiltonian.shape!=self.moving_hamiltonian.shape or self.similarity_hamiltonian.shape!=self.moving_hamiltonian.shape or self.transformation.shape!=self.moving_hamiltonian.shape:
            raise ValueError('lattice oracle matrix shapes are incompatible.')
        for v in (self.maximum_link_unitarity_residual,self.hamiltonian_hermiticity_residual,self.similarity_residual):
            if not np.isfinite(float(v)) or float(v)>tolerance: raise ValueError('lattice moving-frame covariance tolerance failed.')
        return self


def build_lattice_gauge_oracle_v280(model, frame, axes):
    model=model.validate(); frame=frame.validate(); axes=tuple(np.asarray(a,float) for a in axes)
    if len(axes)!=model.ndim or frame.ndim!=model.ndim or frame.nstate!=model.nstate: raise ValueError('lattice oracle dimensions are incompatible.')
    points=lattice_points_v280(axes); kinetic=scalar_lattice_kinetic_v280(axes,model.mass_matrix_au)
    npoint=len(points); nstate=model.nstate; size=npoint*nstate
    fixed=np.kron(kinetic,np.eye(nstate,dtype=complex)); potentials=model.hamiltonian(points)
    for i in range(npoint):
        sl=slice(i*nstate,(i+1)*nstate); fixed[sl,sl]+=potentials[i]
    gauges=frame.unitary(points); transformation=np.zeros((size,size),complex)
    for i in range(npoint):
        sl=slice(i*nstate,(i+1)*nstate); transformation[sl,sl]=gauges[i].conj().T
    moving=np.zeros((size,size),complex); maxlink=0.0
    for i in range(npoint):
        sli=slice(i*nstate,(i+1)*nstate); moving[sli,sli]+=moving_frame_hamiltonian_v280(model,frame,points[i])
        for j in range(npoint):
            value=kinetic[i,j]
            if value==0: continue
            slj=slice(j*nstate,(j+1)*nstate); link=gauges[i].conj().T@gauges[j]
            maxlink=max(maxlink,_scaled_norm_v280(link.conj().T@link,np.eye(nstate))); moving[sli,slj]+=value*link
    similarity=transformation@fixed@transformation.conj().T
    herm=max(_scaled_norm_v280(fixed,fixed.conj().T),_scaled_norm_v280(moving,moving.conj().T)); residual=_scaled_norm_v280(moving,similarity)
    return LatticeGaugeOracleV280(points,fixed,moving,similarity,transformation,maxlink,herm,residual).validate()


def lattice_action_covariance_v280(oracle, reference_state):
    oracle=oracle.validate(); ref=np.asarray(reference_state,complex)
    if ref.shape!=(oracle.fixed_hamiltonian.shape[0],): raise ValueError('lattice state has an incompatible shape.')
    mov=oracle.transformation@ref
    return _scaled_norm_v280(oracle.moving_hamiltonian@mov,oracle.transformation@(oracle.fixed_hamiltonian@ref))


def lattice_propagation_covariance_v280(oracle, reference_state, dt_au):
    oracle=oracle.validate(); ref=np.asarray(reference_state,complex); dt=float(dt_au)
    if ref.shape!=(oracle.fixed_hamiltonian.shape[0],) or not np.isfinite(dt): raise ValueError('lattice propagation inputs are invalid.')
    mov=oracle.transformation@ref; f_end=expm(-1j*dt*oracle.fixed_hamiltonian)@ref; m_end=expm(-1j*dt*oracle.moving_hamiltonian)@mov
    return _scaled_norm_v280(m_end,oracle.transformation@f_end)


def finite_difference_connection_residual_v280(frame, point, step=1e-6):
    frame=frame.validate(); point=np.asarray(point,float); step=float(step)
    if point.shape!=(frame.ndim,) or not np.isfinite(step) or step<=0: raise ValueError('connection finite-difference inputs are invalid.')
    g0=frame.unitary(point); analytic=frame.connection(point); maximum=0.0
    for axis in range(frame.ndim):
        shift=np.zeros(frame.ndim); shift[axis]=step
        deriv=(frame.unitary(point+shift)-frame.unitary(point-shift))/(2*step); numerical=g0.conj().T@deriv
        maximum=max(maximum,_scaled_norm_v280(numerical,analytic[axis]))
    return float(maximum)


def finite_difference_curvature_residual_v280(frame, point, step=2e-5):
    frame=frame.validate(); point=np.asarray(point,float); step=float(step)
    if point.shape!=(frame.ndim,) or frame.ndim<2 or not np.isfinite(step) or step<=0: raise ValueError('curvature finite-difference inputs are invalid.')
    maximum=0.0
    for a,b in itertools.combinations(range(frame.ndim),2):
        ea=np.zeros(frame.ndim); ea[a]=step; eb=np.zeros(frame.ndim); eb[b]=step
        d_a_D_b=(frame.connection(point+ea)[b]-frame.connection(point-ea)[b])/(2*step)
        d_b_D_a=(frame.connection(point+eb)[a]-frame.connection(point-eb)[a])/(2*step)
        D=frame.connection(point); curvature=d_a_D_b-d_b_D_a+D[a]@D[b]-D[b]@D[a]
        maximum=max(maximum,float(np.linalg.norm(curvature,ord='fro')))
    return maximum

__all__=['LatticeGaugeOracleV280','build_lattice_gauge_oracle_v280','finite_difference_connection_residual_v280','finite_difference_curvature_residual_v280','lattice_action_covariance_v280','lattice_points_v280','lattice_propagation_covariance_v280','periodic_second_derivative_v280','scalar_lattice_kinetic_v280']
