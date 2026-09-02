from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


@dataclass
class SparseSpinorMatrices:
    S: sparse.csr_matrix
    H: sparse.csr_matrix
    Snuc: sparse.csr_matrix
    active_edges: tuple
    n_basis: int
    nstate: int = 2

    @property
    def electronic_dimension(self):
        return self.nstate*self.n_basis

    @property
    def S_density(self):
        n=self.S.shape[0]
        return float(self.S.nnz/max(n*n,1))

    @property
    def H_density(self):
        n=self.H.shape[0]
        return float(self.H.nnz/max(n*n,1))

    @property
    def nuclear_density(self):
        n=self.Snuc.shape[0]
        return float(self.Snuc.nnz/max(n*n,1))

    def as_dict(self):
        return {
            "n_basis":int(self.n_basis),
            "electronic_dimension":
                int(self.electronic_dimension),
            "active_offdiagonal_edges":
                int(len(self.active_edges)),
            "S_nnz":int(self.S.nnz),
            "H_nnz":int(self.H.nnz),
            "Snuc_nnz":int(self.Snuc.nnz),
            "S_density":self.S_density,
            "H_density":self.H_density,
            "Snuc_density":self.nuclear_density,
        }


def _append_block(rows,cols,data,i,j,block,nstate=2):
    block=np.asarray(block,dtype=complex)
    if block.shape!=(nstate,nstate):
        raise ValueError("block has incompatible shape.")
    for a in range(nstate):
        for b in range(nstate):
            value=block[a,b]
            if value!=0.0:
                rows.append(nstate*i+a)
                cols.append(nstate*j+b)
                data.append(value)


def build_sparse_spinor_lvc_matrices(
    locality_update,
    provider,
    *,
    nstate=2,
):
    r"""Build sparse exact-LVC S/H matrices on the active locality graph.

    Diagonal Gaussian blocks are always retained.

    For an active undirected edge (i,j), the exact Gaussian pair data is evaluated and
    both Hermitian orientations are inserted.

    Pairs absent from the locality graph are represented by exact structural zeros in
    the sparse approximation.
    """
    if int(nstate)!=2:
        raise ValueError("analytic LVC sparse builder currently requires nstate=2.")

    cache=locality_update.cache
    basis=cache.basis
    n=len(basis)
    dim=nstate*n

    point=provider.evaluate(
        np.asarray(basis[0].q,float)
    )
    Minv=np.linalg.inv(
        np.asarray(point.mass_matrix,float)
    )
    params=provider.params
    eye=np.eye(nstate,dtype=complex)

    s_rows=[]; s_cols=[]; s_data=[]
    h_rows=[]; h_cols=[]; h_data=[]
    n_rows=[]; n_cols=[]; n_data=[]

    # Diagonal blocks.
    for i in range(n):
        pair=cache.pair(i,i)
        Sii=pair.overlap
        Kii=pair.kinetic(Minv)
        Vii=pair.lvc_potential_matrix(params)

        _append_block(
            s_rows,s_cols,s_data,
            i,i,Sii*eye,nstate
        )
        _append_block(
            h_rows,h_cols,h_data,
            i,i,Kii*eye+Vii,nstate
        )
        n_rows.append(i); n_cols.append(i); n_data.append(Sii)

    # Active off-diagonal graph.
    for i,j in locality_update.active_edges:
        pair=cache.pair(i,j)
        Sij=pair.overlap
        Kij=pair.kinetic(Minv)
        Vij=pair.lvc_potential_matrix(params)

        blockS=Sij*eye
        blockH=Kij*eye+Vij

        _append_block(
            s_rows,s_cols,s_data,
            i,j,blockS,nstate
        )
        _append_block(
            s_rows,s_cols,s_data,
            j,i,blockS.conj().T,nstate
        )
        _append_block(
            h_rows,h_cols,h_data,
            i,j,blockH,nstate
        )
        _append_block(
            h_rows,h_cols,h_data,
            j,i,blockH.conj().T,nstate
        )

        n_rows.extend([i,j])
        n_cols.extend([j,i])
        n_data.extend([Sij,np.conj(Sij)])

    S=sparse.coo_matrix(
        (np.asarray(s_data,dtype=complex),(s_rows,s_cols)),
        shape=(dim,dim),
    ).tocsr()
    H=sparse.coo_matrix(
        (np.asarray(h_data,dtype=complex),(h_rows,h_cols)),
        shape=(dim,dim),
    ).tocsr()
    Snuc=sparse.coo_matrix(
        (np.asarray(n_data,dtype=complex),(n_rows,n_cols)),
        shape=(n,n),
    ).tocsr()

    return SparseSpinorMatrices(
        S=S,H=H,Snuc=Snuc,
        active_edges=tuple(locality_update.active_edges),
        n_basis=n,
        nstate=nstate,
    )


def build_sparse_spinor_time_matrix(
    locality_update,
    qdots,
    pdots,
    Adots=None,
    *,
    nstate=2,
):
    r"""Build sparse ordered moving-basis T only on diagonal + active graph edges."""
    if int(nstate)!=2:
        raise ValueError("current spinor-complete sparse builder requires nstate=2.")

    cache=locality_update.cache
    n=len(cache)
    dim=nstate*n
    nq=cache.dimension

    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)
    if qdots.shape!=(n,nq) or pdots.shape!=qdots.shape:
        raise ValueError("qdots/pdots have incompatible shape.")

    if Adots is None:
        Adots=[None]*n
    if len(Adots)!=n:
        raise ValueError("Adots must contain one entry per Gaussian.")

    rows=[]; cols=[]; data=[]
    eye=np.eye(nstate,dtype=complex)

    def add_oriented(i,j):
        pair=cache.pair(i,j)
        tij=pair.time_element(
            qdots[j],pdots[j],Adots[j]
        )
        _append_block(
            rows,cols,data,
            i,j,tij*eye,nstate
        )

    for i in range(n):
        add_oriented(i,i)

    for i,j in locality_update.active_edges:
        add_oriented(i,j)
        add_oriented(j,i)

    return sparse.coo_matrix(
        (np.asarray(data,dtype=complex),(rows,cols)),
        shape=(dim,dim),
    ).tocsr()


def sparse_metric_compatible_connection(
    S_old,
    S_new,
    dt,
    seed,
):
    r"""Sparse version of T = T0 + 1/2[dS - T0 - T0^dagger]."""
    dt=float(dt)
    if dt<=0.0:
        raise ValueError("dt must be positive.")

    S0=sparse.csr_matrix(S_old,dtype=complex)
    S1=sparse.csr_matrix(S_new,dtype=complex)
    T0=sparse.csr_matrix(seed,dtype=complex)

    if S0.shape!=S1.shape or T0.shape!=S0.shape:
        raise ValueError("S_old/S_new/seed shapes must match.")

    dS=(S1-S0)*(1.0/dt)
    correction=0.5*(
        dS-T0-T0.getH()
    )
    return (T0+correction).tocsr()


def sparse_moving_basis_midpoint_cayley_step(
    C,
    S_old,
    H_old,
    S_new,
    H_new,
    T_mid,
    dt,
):
    r"""Sparse implicit midpoint/Cayley solve for i S Cdot=(H-iT)C."""
    C=np.asarray(C,dtype=complex)
    dt=float(dt)
    if dt<=0.0:
        raise ValueError("dt must be positive.")

    S0=sparse.csr_matrix(S_old,dtype=complex)
    S1=sparse.csr_matrix(S_new,dtype=complex)
    H0=sparse.csr_matrix(H_old,dtype=complex)
    H1=sparse.csr_matrix(H_new,dtype=complex)
    T=sparse.csr_matrix(T_mid,dtype=complex)

    if not (
        S0.shape==S1.shape==H0.shape==H1.shape==T.shape
    ):
        raise ValueError("all sparse matrices must have equal shape.")
    if C.shape!=(S0.shape[0],):
        raise ValueError("coefficient vector has incompatible size.")

    Sm=0.5*(S0+S1)
    Hm=0.5*(H0+H1)
    K=1j*Hm+T

    lhs=(Sm+0.5*dt*K).tocsc()
    rhs=(Sm-0.5*dt*K)@C
    return np.asarray(
        spsolve(lhs,rhs),
        dtype=complex,
    )


def sparse_generalized_norm(C,S):
    C=np.asarray(C,dtype=complex)
    S=sparse.csr_matrix(S,dtype=complex)
    return float(np.real(np.vdot(C,S@C)))


def sparse_reduced_density(flat_coefficients,Snuc,nstate=2,normalize=True):
    S=sparse.csr_matrix(Snuc,dtype=complex)
    n=S.shape[0]
    C=np.asarray(flat_coefficients,dtype=complex).reshape(n,nstate)

    # rho_ab = sum_ij C_ia S_ij C_jb^*
    SC=S@C
    rho=C.T@np.conj(SC)
    rho=0.5*(rho+rho.conj().T)

    if normalize:
        tr=np.trace(rho)
        if abs(tr)<1e-15:
            raise ValueError("zero reduced-density trace.")
        rho=rho/tr
    return rho


def sparse_matrix_relative_difference(A,B,ord="fro"):
    """Reference utility for small validation problems."""
    Ad=sparse.csr_matrix(A).toarray()
    Bd=sparse.csr_matrix(B).toarray()
    denom=max(np.linalg.norm(Bd,ord=ord),1e-30)
    return float(np.linalg.norm(Ad-Bd,ord=ord)/denom)


def audit_sparse_lvc_matrices_against_dense(
    basis,
    provider,
    sparse_mats,
):
    r"""A-posteriori dense audit of one sparse snapshot.

    This intentionally pays the O(N^2 d^3) dense reference cost at an audit checkpoint.
    It is not part of every sparse propagation step.

    The audit quantifies the approximation introduced by dropping locality-graph
    edges.  In particular, an overlap cutoff by itself is not a rigorous upper bound
    on kinetic/potential matrix-element error because polynomial prefactors can
    amplify a small overlap.  The dense audit therefore reports S and H errors
    explicitly.
    """
    from .pair_cache_v15 import (
        GaussianPairCache,
        build_cached_spinor_lvc_matrices,
    )

    cache=GaussianPairCache(basis)
    Sd,Hd,Nd=build_cached_spinor_lvc_matrices(
        cache,provider
    )

    Ss=sparse_mats.S.toarray()
    Hs=sparse_mats.H.toarray()
    Ns=sparse_mats.Snuc.toarray()

    def rel(a,b):
        return float(
            np.linalg.norm(a-b,ord="fro")
            /max(np.linalg.norm(b,ord="fro"),1e-30)
        )

    omitted_S=Sd-Ss
    omitted_H=Hd-Hs

    n=len(basis)
    max_omitted_overlap=0.0
    max_omitted_h_block=0.0
    omitted_edge_count=0

    active=set(tuple(x) for x in sparse_mats.active_edges)
    for i in range(n):
        for j in range(i+1,n):
            if (i,j) in active:
                continue
            omitted_edge_count+=1
            max_omitted_overlap=max(
                max_omitted_overlap,
                float(abs(Nd[i,j])),
            )
            si=slice(2*i,2*i+2)
            sj=slice(2*j,2*j+2)
            max_omitted_h_block=max(
                max_omitted_h_block,
                float(np.linalg.norm(
                    Hd[si,sj],ord="fro"
                )),
            )

    return {
        "relative_S_frobenius_error":rel(Ss,Sd),
        "relative_H_frobenius_error":rel(Hs,Hd),
        "relative_Snuc_frobenius_error":rel(Ns,Nd),
        "maximum_omitted_overlap":
            float(max_omitted_overlap),
        "maximum_omitted_H_block_frobenius":
            float(max_omitted_h_block),
        "omitted_offdiagonal_edges":
            int(omitted_edge_count),
        "dense_pair_factorizations":
            int(cache.stats.canonical_solves),
    }
