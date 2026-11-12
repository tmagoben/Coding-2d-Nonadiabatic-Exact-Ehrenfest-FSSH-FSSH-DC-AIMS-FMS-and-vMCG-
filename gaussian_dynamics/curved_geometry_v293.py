r"""Coordinate-dependent non-Abelian geometry for v0.29.3 development.

The connection convention is D_mu = partial_mu + A_mu with A_mu anti-Hermitian.
Wilson transport W(b<-a;gamma) solves

    dW/ds = - dot(R)^mu A_mu(R(s)) W,   W(0)=I.

Under a coordinate-dependent unitary frame change G(R),

    A'_mu = G^\dagger A_mu G + G^\dagger partial_mu G,
    W'(b<-a) = G(b)^\dagger W(b<-a) G(a).

This module deliberately treats the path as part of the section definition.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

Array = np.ndarray


def _unitary_error(u: Array) -> float:
    u=np.asarray(u,dtype=complex)
    return float(np.linalg.norm(u.conj().T@u-np.eye(u.shape[0])))


@dataclass(frozen=True)
class PolylinePathV293:
    vertices: Array
    label: str = "polyline"
    def __post_init__(self):
        object.__setattr__(self,"vertices",np.asarray(self.vertices,dtype=float).copy())
    @property
    def ndim(self): return int(self.vertices.shape[1]) if self.vertices.ndim==2 else 0
    @property
    def start(self): return self.vertices[0].copy()
    @property
    def end(self): return self.vertices[-1].copy()
    def validate(self):
        v=self.vertices
        if v.ndim!=2 or v.shape[0]<2 or v.shape[1]<1 or not np.all(np.isfinite(v)):
            raise ValueError("path vertices must have shape (n>=2,ndim) and be finite")
        constant_path = len(v) == 2 and np.array_equal(v[0], v[1])
        if not constant_path and any(np.linalg.norm(v[i+1]-v[i]) <= 1e-15 for i in range(len(v)-1)):
            raise ValueError("path contains a zero-length segment")
        return self
    def reversed(self):
        return PolylinePathV293(self.vertices[::-1],label=f"reverse({self.label})").validate()


def straight_path_v293(start, end, *, label="straight"):
    return PolylinePathV293(np.vstack([np.asarray(start,float),np.asarray(end,float)]),label).validate()


def dogleg_path_v293(start,end,*,axis=0,fraction=0.5,label="dogleg"):
    a=np.asarray(start,float); b=np.asarray(end,float)
    if a.shape!=b.shape or a.ndim!=1: raise ValueError("path endpoints disagree")
    axis=int(axis)
    if not (0<=axis<a.size): raise ValueError("invalid dogleg axis")
    mid=a.copy(); mid[axis]=a[axis]+float(fraction)*(b[axis]-a[axis])
    # then move every other coordinate, then finish chosen axis
    mid2=b.copy(); mid2[axis]=mid[axis]
    verts=[a]
    if np.linalg.norm(mid-a)>1e-15: verts.append(mid)
    if np.linalg.norm(mid2-verts[-1])>1e-15: verts.append(mid2)
    if np.linalg.norm(b-verts[-1])>1e-15: verts.append(b)
    if len(verts)<2: verts=[a,b]
    return PolylinePathV293(np.asarray(verts),label).validate()


@dataclass(frozen=True)
class ConnectionFieldV293:
    ndim: int
    nstate: int
    evaluator: Callable[[Array], Array]
    name: str = "connection-field"
    segment_transport: Callable[[Array,Array],Array] | None = None
    def at(self, point):
        r=np.asarray(point,dtype=float)
        if r.shape!=(self.ndim,) or not np.all(np.isfinite(r)): raise ValueError("connection point has wrong shape")
        a=np.asarray(self.evaluator(r),dtype=complex)
        if a.shape!=(self.ndim,self.nstate,self.nstate) or not np.all(np.isfinite(a)):
            raise ValueError("connection evaluator returned wrong shape/non-finite values")
        for m in range(self.ndim):
            scale=max(float(np.linalg.norm(a[m])),1.0)
            if np.linalg.norm(a[m]+a[m].conj().T)/scale>5e-10:
                raise ValueError("connection must be anti-Hermitian")
        return a
    def validate(self, probes=None):
        if self.ndim<1 or self.nstate<1: raise ValueError("invalid connection dimensions")
        if probes is not None:
            for r in probes: self.at(r)
        return self
    def curvature(self, point, mu, nu, *, step=2e-5):
        r=np.asarray(point,float); mu=int(mu); nu=int(nu)
        if mu==nu: return np.zeros((self.nstate,self.nstate),complex)
        em=np.zeros(self.ndim); en=np.zeros(self.ndim); em[mu]=step; en[nu]=step
        da_nu=(self.at(r+em)[nu]-self.at(r-em)[nu])/(2*step)
        da_mu=(self.at(r+en)[mu]-self.at(r-en)[mu])/(2*step)
        a=self.at(r)
        return da_nu-da_mu+a[mu]@a[nu]-a[nu]@a[mu]


def constant_connection_field_v293(generators, *, name="constant-connection"):
    g=np.asarray(generators,dtype=complex)
    if g.ndim!=3 or g.shape[1]!=g.shape[2]: raise ValueError("generators must have shape (ndim,nstate,nstate)")
    nd,ns=g.shape[:2]
    def eval_(_r): return g
    def transport_(a,b):
        dr=np.asarray(b,float)-np.asarray(a,float)
        return expm(-np.einsum("m,mab->ab",dr,g,optimize=True))
    return ConnectionFieldV293(nd,ns,eval_,name,transport_).validate([np.zeros(nd)])


@dataclass(frozen=True)
class LinearConnectionFieldV293:
    base: Array
    gradients: Array
    name: str = "linear-connection"
    def __post_init__(self):
        object.__setattr__(self,"base",np.asarray(self.base,dtype=complex).copy())
        object.__setattr__(self,"gradients",np.asarray(self.gradients,dtype=complex).copy())
    @property
    def ndim(self): return int(self.base.shape[0]) if self.base.ndim==3 else 0
    @property
    def nstate(self): return int(self.base.shape[1]) if self.base.ndim==3 else 0
    def field(self):
        base=self.base.copy(); grad=self.gradients.copy(); nd=self.ndim; ns=self.nstate
        if grad.shape!=(nd,nd,ns,ns): raise ValueError("gradients must have shape (coordinate,component,state,state)")
        def eval_(r): return base+np.einsum('x,xmab->mab',r,grad,optimize=True)
        return ConnectionFieldV293(nd,ns,eval_,self.name).validate([np.zeros(nd)])


@dataclass(frozen=True)
class CoordinateGaugeV293:
    ndim: int
    nstate: int
    matrix: Callable[[Array],Array]
    derivative: Callable[[Array],Array] | None = None
    name: str = "gauge"
    def at(self,r):
        u=np.asarray(self.matrix(np.asarray(r,float)),complex)
        if u.shape!=(self.nstate,self.nstate) or not np.all(np.isfinite(u)) or _unitary_error(u)>5e-9: raise ValueError("gauge must be finite and unitary")
        return u
    def deriv(self,r,*,step=2e-6):
        r=np.asarray(r,float)
        if self.derivative is not None:
            d=np.asarray(self.derivative(r),complex)
            if d.shape!=(self.ndim,self.nstate,self.nstate) or not np.all(np.isfinite(d)): raise ValueError("gauge derivative has wrong shape or non-finite values")
            return d
        out=[]
        for mu in range(self.ndim):
            e=np.zeros(self.ndim); e[mu]=step
            out.append((self.at(r+e)-self.at(r-e))/(2*step))
        return np.asarray(out)


def gauge_transform_connection_v293(connection: ConnectionFieldV293, gauge: CoordinateGaugeV293):
    connection.validate();
    if (connection.ndim,connection.nstate)!=(gauge.ndim,gauge.nstate): raise ValueError("gauge/connection dimensions differ")
    def eval_(r):
        g=gauge.at(r); dg=gauge.deriv(r); a=connection.at(r)
        return np.asarray([g.conj().T@a[mu]@g + g.conj().T@dg[mu] for mu in range(connection.ndim)])
    def transport_(a,b):
        if connection.segment_transport is not None:
            w=connection.segment_transport(np.asarray(a,float),np.asarray(b,float))
        else:
            w=wilson_path_v293(connection,straight_path_v293(a,b))
        return gauge.at(b).conj().T @ w @ gauge.at(a)
    return ConnectionFieldV293(connection.ndim,connection.nstate,eval_,f"{connection.name}|{gauge.name}",transport_)


def wilson_path_v293(connection: ConnectionFieldV293, path: PolylinePathV293, *, rtol=2e-11, atol=2e-13):
    connection.validate(); path=path.validate()
    if path.ndim!=connection.ndim: raise ValueError("path/connection dimensions differ")
    w=np.eye(connection.nstate,dtype=complex)
    for a,b in zip(path.vertices[:-1],path.vertices[1:]):
        if np.array_equal(a, b):
            continue
        if connection.segment_transport is not None:
            segment=np.asarray(connection.segment_transport(a,b),complex)
            if segment.shape!=(connection.nstate,connection.nstate) or not np.all(np.isfinite(segment)):
                raise ValueError("segment transport has wrong shape or non-finite values")
            if _unitary_error(segment)>5e-9:
                raise ValueError("segment transport must be unitary")
            w=segment @ w
            continue
        dr=b-a
        def rhs(s,y):
            mat=y.reshape(connection.nstate,connection.nstate)
            r=a+s*dr
            generator=-np.einsum('m,mab->ab',dr,connection.at(r),optimize=True)
            return (generator@mat).reshape(-1)
        sol=solve_ivp(rhs,(0.0,1.0),w.reshape(-1),method="DOP853",rtol=rtol,atol=atol)
        if not sol.success: raise RuntimeError("Wilson ODE failed: "+sol.message)
        w=sol.y[:,-1].reshape(connection.nstate,connection.nstate)
    return w


def wilson_straight_v293(connection,start,end,**kwargs):
    return wilson_path_v293(connection,straight_path_v293(start,end),**kwargs)


def plaquette_holonomy_v293(connection, origin, mu, nu, h):
    r=np.asarray(origin,float); em=np.zeros(connection.ndim); en=np.zeros(connection.ndim); em[mu]=h; en[nu]=h
    p=PolylinePathV293(np.vstack([r,r+em,r+em+en,r+en,r]),label="plaquette").validate()
    return wilson_path_v293(connection,p)


@dataclass(frozen=True)
class PatchTransitionV293:
    source_name: str
    target_name: str
    gauge: CoordinateGaugeV293  # target frame = source frame * G
    def matrix_target_from_source(self,r):
        return self.gauge.at(r).conj().T


def validate_patch_transition_v293(source: ConnectionFieldV293,target: ConnectionFieldV293,transition: PatchTransitionV293,probes,*,tol=2e-7):
    transformed=gauge_transform_connection_v293(source,transition.gauge)
    worst=0.0
    for r in probes:
        a=transformed.at(r); b=target.at(r)
        worst=max(worst,float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1.0)))
    if worst>tol: raise ValueError("patch transition does not relate the two local connections")
    return worst


def two_patch_transport_v293(source,target,transition,start,switch,end,*,source_path_factory=straight_path_v293,target_path_factory=straight_path_v293):
    ws=wilson_path_v293(source,source_path_factory(start,switch))
    wt=wilson_path_v293(target,target_path_factory(switch,end))
    return wt @ transition.matrix_target_from_source(switch) @ ws


def section_transport_v293(connection,start,end,path_factory=straight_path_v293):
    return wilson_path_v293(connection,path_factory(start,end))
