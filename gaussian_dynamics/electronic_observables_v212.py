from dataclasses import dataclass
import numpy as np
from scipy import sparse

from .gauge_graph import nearest_unitary
from .gaussian_general import gaussian_overlap_general, real_overlap_saddle_point
from .matrix_invariants_v213 import hermiticity_residual_v213, require_residual_v213


@dataclass(frozen=True)
class ElectronicObservableV212:
    """Representation-neutral electronic observable evaluated in a local frame.

    evaluator(snapshot) must return the matrix of one physical Hermitian electronic
    operator in the electronic frame of `snapshot`.  Examples can later include spin,
    diabatic character, charge, or any other electronic operator without changing the
    expectation-value algebra.
    """
    name: str
    evaluator: object

    def matrix(self,snapshot):
        O=np.asarray(self.evaluator(snapshot),dtype=complex)
        ns=snapshot.point.nstate
        if O.shape!=(ns,ns):
            raise ValueError("observable matrix has incompatible shape.")
        require_residual_v213(
            f"observable '{self.name}' Hermiticity",
            hermiticity_residual_v213(O),
            1.0e-12,
        )
        return O


def build_electronic_observable_matrix_v212(basis, provider, observable, active_edges=None):
    basis=list(basis)
    if not basis:
        raise ValueError("basis cannot be empty.")
    centers=[provider.evaluate_snapshot(b.q) for b in basis]
    s=centers[0].point.nstate
    n=len(basis)
    dim=n*s
    Omat=np.zeros((dim,dim),dtype=complex)

    def put(i,j,B):
        Omat[s*i:s*(i+1),s*j:s*(j+1)]=B

    for i,snap in enumerate(centers):
        put(i,i,observable.matrix(snap))

    if active_edges is None:
        edges=((i,j) for i in range(n) for j in range(i+1,n))
    else:
        edges=tuple(tuple(map(int,e)) for e in active_edges)

    for i,j in edges:
        bi,bj=basis[i],basis[j]
        qbar=real_overlap_saddle_point(bi.q,bi.A,bj.q,bj.A)
        sc=provider.evaluate_snapshot(qbar)
        Uci=nearest_unitary(np.asarray(provider.snapshot_overlap(sc,centers[i]),complex))
        Ucj=nearest_unitary(np.asarray(provider.snapshot_overlap(sc,centers[j]),complex))
        Oc=observable.matrix(sc)
        block_e=Uci.conj().T@Oc@Ucj
        Sn=gaussian_overlap_general(bi.q,bi.p,bi.A,bj.q,bj.p,bj.A)
        block=Sn*block_e
        put(i,j,block)
        put(j,i,block.conj().T)

    return sparse.csr_matrix(Omat)


def observable_expectation_v212(coefficients, metric, observable_matrix):
    C=np.asarray(coefficients,dtype=complex)
    S=metric
    O=observable_matrix
    denom=np.vdot(C,S@C)
    if abs(denom)<1e-30:
        raise ValueError("represented state has zero metric norm.")
    value=np.vdot(C,O@C)/denom
    if abs(np.imag(value))>1e-8:
        raise ValueError("Hermitian observable expectation acquired a significant imaginary part.")
    return float(np.real(value))
