from dataclasses import dataclass, field, replace
import math
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from .gaussian_nd import analytic_overlap_equal_width
from .gaussian_general import validate_spd
from .gauge_graph import nearest_unitary
from .local_gaussian_nd import (
    kinetic_matrix_element_equal_width,
    basis_time_matrix_element_equal_width,
)


@dataclass
class SparseMolecularTBFV20:
    uid: int
    state: int
    q: np.ndarray
    p: np.ndarray
    A: np.ndarray

    def __post_init__(self):
        self.uid=int(self.uid)
        self.state=int(self.state)
        self.q=np.asarray(self.q,dtype=float)
        self.p=np.asarray(self.p,dtype=float)
        self.A=validate_spd(self.A)
        if self.q.shape!=self.p.shape:
            raise ValueError("q and p must have equal shapes.")
        if self.A.shape!=(len(self.q),len(self.q)):
            raise ValueError("A has incompatible shape.")

    def copy(self):
        return SparseMolecularTBFV20(
            self.uid,self.state,
            self.q.copy(),self.p.copy(),self.A.copy(),
        )


@dataclass(frozen=True)
class MolecularSparseSettingsV20:
    enter_score: float=0.030
    exit_score: float=0.015
    search_overlap_floor: float=1e-5
    overlap_weight: float=1.0
    hamiltonian_weight: float=0.20
    time_connection_weight: float=1.0
    local_omitted_score_l2_budget: float=0.010
    hamiltonian_floor: float=1e-10
    use_kdtree: bool=True

    def validate(self):
        if not (0.0<self.exit_score<=self.enter_score):
            raise ValueError("Require 0 < exit_score <= enter_score.")
        if not (0.0<self.search_overlap_floor<1.0):
            raise ValueError("search_overlap_floor must lie in (0,1).")
        if self.local_omitted_score_l2_budget<0.0:
            raise ValueError("local omitted-score budget cannot be negative.")
        if self.hamiltonian_floor<=0.0:
            raise ValueError("hamiltonian_floor must be positive.")
        return self


@dataclass(frozen=True)
class MolecularPairDataV20:
    i: int
    j: int
    nuclear_overlap: complex
    electronic_overlap: complex
    S: complex
    H: complex
    T_ij: complex
    T_ji: complex
    hamiltonian_relative: float
    time_connection_dt: float
    score: float
    centroid_q: np.ndarray


@dataclass
class MolecularSparseUpdateV20:
    active_edges: tuple
    pair_data: dict
    diagonal_S: np.ndarray
    diagonal_H: np.ndarray
    diagonal_T: np.ndarray
    qdots: np.ndarray
    pdots: np.ndarray
    spatial_candidate_pairs: int
    candidate_pairs: int
    exact_pair_checks: int
    budget_promoted_edges: int
    omitted_score_l2: float
    total_offdiagonal_pairs: int
    entered_edges: int
    exited_edges: int
    retained_edges: int
    maximum_omitted_candidate_score: float

    @property
    def active_offdiagonal_edges(self):
        return len(self.active_edges)

    @property
    def sparsity_fraction(self):
        total=max(self.total_offdiagonal_pairs,1)
        return float(
            1.0-len(self.active_edges)/total
        )

    def as_dict(self):
        return {
            "active_offdiagonal_edges":
                int(self.active_offdiagonal_edges),
            "spatial_candidate_pairs":
                int(self.spatial_candidate_pairs),
            "candidate_pairs":
                int(self.candidate_pairs),
            "exact_pair_checks":
                int(self.exact_pair_checks),
            "budget_promoted_edges":
                int(self.budget_promoted_edges),
            "omitted_score_l2":
                float(self.omitted_score_l2),
            "total_offdiagonal_pairs":
                int(self.total_offdiagonal_pairs),
            "entered_edges":
                int(self.entered_edges),
            "exited_edges":
                int(self.exited_edges),
            "retained_edges":
                int(self.retained_edges),
            "maximum_omitted_candidate_score":
                float(self.maximum_omitted_candidate_score),
            "sparsity_fraction":
                self.sparsity_fraction,
        }


@dataclass
class SparseMolecularMatricesV20:
    S: sparse.csr_matrix
    H: sparse.csr_matrix
    T_seed: sparse.csr_matrix
    active_edges: tuple
    update: MolecularSparseUpdateV20


def _uid_edge(a,b):
    a=int(a); b=int(b)
    return (a,b) if a<b else (b,a)


def _position_bound(qi,qj,ai,aj):
    h=1.0/(1.0/float(ai)+1.0/float(aj))
    dq=np.asarray(qi,float)-np.asarray(qj,float)
    return float(
        math.exp(-0.5*h*float(dq@dq))
    )


def _safe_radius(amin,threshold):
    return float(math.sqrt(
        max(
            -4.0*math.log(float(threshold))
            /float(amin),
            0.0,
        )
    ))


def _state_unit(nstate,state):
    c=np.zeros(int(nstate),dtype=complex)
    c[int(state)]=1.0
    return c


def _center_kinematics(basis,provider):
    qdots=[]; pdots=[]; snapshots=[]
    for b in basis:
        snap=provider.evaluate_snapshot(b.q)
        point=snap.point
        qdot=np.linalg.solve(
            point.mass_matrix_q_au,b.p
        )
        pdot=-point.gradients_q[b.state]
        qdots.append(qdot)
        pdots.append(pdot)
        snapshots.append(snap)
    return (
        np.asarray(qdots,float),
        np.asarray(pdots,float),
        snapshots,
    )


def _transport_to_centroid(
    provider,
    center_snapshot,
    centroid_snapshot,
    state,
):
    O=provider.snapshot_overlap(
        centroid_snapshot,
        center_snapshot,
    )
    U=nearest_unitary(O)
    return U@_state_unit(
        len(centroid_snapshot.point.energies),
        state,
    )


def molecular_pair_data_v20(
    basis,
    i,
    j,
    provider,
    dt,
    qdots,
    pdots,
    center_snapshots,
    diagonal_h_abs,
    settings,
):
    """Exact active-edge pair data for the v0.20 centroid approximation."""
    i=int(i); j=int(j)
    if i==j:
        raise ValueError("molecular_pair_data_v20 is for off-diagonal pairs.")

    bi=basis[i]; bj=basis[j]
    if not np.allclose(
        bi.A,bj.A,atol=1e-12
    ):
        raise ValueError(
            "v0.20 molecular pair approximation currently requires equal widths."
        )

    qbar=0.5*(bi.q+bj.q)
    centroid=provider.evaluate_snapshot(qbar)
    M=centroid.point.mass_matrix_q_au

    vi=_transport_to_centroid(
        provider,center_snapshots[i],
        centroid,bi.state,
    )
    vj=_transport_to_centroid(
        provider,center_snapshots[j],
        centroid,bj.state,
    )

    e_overlap=np.vdot(vi,vj)
    He=np.diag(
        np.asarray(
            centroid.point.energies,
            dtype=float,
        )
    ).astype(complex)
    potential=np.vdot(
        vi,He@vj
    )

    Snuc=analytic_overlap_equal_width(
        bi.q,bi.p,bj.q,bj.p,bi.A
    )
    Tnuc=kinetic_matrix_element_equal_width(
        bi.q,bi.p,bj.q,bj.p,bi.A,M
    )

    S=Snuc*e_overlap
    H=Tnuc*e_overlap+Snuc*potential

    Tn_ij=basis_time_matrix_element_equal_width(
        bi.q,bi.p,bj.q,bj.p,bj.A,
        qdots[j],pdots[j],
    )
    Tn_ji=basis_time_matrix_element_equal_width(
        bj.q,bj.p,bi.q,bi.p,bi.A,
        qdots[i],pdots[i],
    )
    T_ij=Tn_ij*e_overlap
    T_ji=Tn_ji*np.conj(e_overlap)

    hscale=max(
        math.sqrt(
            max(float(diagonal_h_abs[i]),0.0)
            *max(float(diagonal_h_abs[j]),0.0)
        ),
        settings.hamiltonian_floor,
    )
    hrel=float(abs(H)/hscale)
    tdt=float(
        float(dt)*math.sqrt(
            abs(T_ij)**2+abs(T_ji)**2
        )
    )
    score=float(math.sqrt(
        (
            settings.overlap_weight*abs(S)
        )**2
        +(
            settings.hamiltonian_weight*hrel
        )**2
        +(
            settings.time_connection_weight*tdt
        )**2
    ))

    return MolecularPairDataV20(
        i=i,j=j,
        nuclear_overlap=Snuc,
        electronic_overlap=e_overlap,
        S=S,H=H,
        T_ij=T_ij,T_ji=T_ji,
        hamiltonian_relative=hrel,
        time_connection_dt=tdt,
        score=score,
        centroid_q=qbar.copy(),
    )


class SparseMolecularEdgeGraphV20:
    """Persistent molecular S/H/T importance graph with geometric pre-screening."""

    def __init__(
        self,
        provider,
        dt,
        settings=MolecularSparseSettingsV20(),
    ):
        self.provider=provider
        self.dt=float(dt)
        if self.dt<=0.0:
            raise ValueError("dt must be positive.")
        self.settings=settings.validate()
        self._active_uid_edges=set()
        self.update_count=0
        self.total_exact_pair_checks=0

    @property
    def active_uid_edges(self):
        return tuple(sorted(self._active_uid_edges))

    def relax_search_floor(self,factor):
        factor=float(factor)
        if not (0.0<factor<1.0):
            raise ValueError("relaxation factor must lie in (0,1).")
        self.settings=replace(
            self.settings,
            search_overlap_floor=max(
                self.settings.search_overlap_floor*factor,
                1e-14,
            ),
        )
        return self.settings

    def update(self,basis):
        basis=list(basis)
        if not basis:
            raise ValueError("basis cannot be empty.")
        uids=[int(b.uid) for b in basis]
        if len(set(uids))!=len(uids):
            raise ValueError("TBF uids must be unique.")
        uid_to_index={
            uid:i for i,uid in enumerate(uids)
        }
        live=set(uids)
        old={
            e for e in self._active_uid_edges
            if e[0] in live and e[1] in live
        }

        q=np.asarray(
            [b.q for b in basis],float
        )
        amin=np.asarray([
            float(np.min(
                np.linalg.eigvalsh(
                    validate_spd(b.A)
                )
            ))
            for b in basis
        ])

        qdots,pdots,centers=(
            _center_kinematics(
                basis,self.provider
            )
        )

        # Diagonal blocks need no pair centroid.
        n=len(basis)
        diag_S=np.ones(n,dtype=complex)
        diag_H=np.zeros(n,dtype=complex)
        diag_T=np.zeros(n,dtype=complex)
        for i,b in enumerate(basis):
            point=centers[i].point
            Tii=kinetic_matrix_element_equal_width(
                b.q,b.p,b.q,b.p,b.A,
                point.mass_matrix_q_au,
            )
            diag_H[i]=(
                Tii+point.energies[b.state]
            )
            diag_T[i]=(
                basis_time_matrix_element_equal_width(
                    b.q,b.p,b.q,b.p,b.A,
                    qdots[i],pdots[i],
                )
            )
        diag_h_abs=np.abs(diag_H)

        total=n*(n-1)//2
        if (
            self.settings.use_kdtree
            and n>1
        ):
            radius=_safe_radius(
                np.min(amin),
                self.settings.search_overlap_floor,
            )
            candidates=set(
                tuple(sorted(x))
                for x in cKDTree(q).query_pairs(
                    radius,
                    output_type="set",
                )
            )
        else:
            candidates={
                (i,j)
                for i in range(n)
                for j in range(i+1,n)
            }

        # Existing edges are always exactly rescored so hysteresis cannot strand them.
        for ua,ub in old:
            candidates.add(
                tuple(sorted((
                    uid_to_index[ua],
                    uid_to_index[ub],
                )))
            )

        spatial_count=len(candidates)

        # Pair-specific conservative nuclear-overlap search floor.
        filtered=[]
        for i,j in sorted(candidates):
            edge=_uid_edge(
                uids[i],uids[j]
            )
            if edge in old:
                filtered.append((i,j))
                continue
            upper=_position_bound(
                q[i],q[j],amin[i],amin[j]
            )
            if (
                upper
                >=self.settings.search_overlap_floor
            ):
                filtered.append((i,j))

        pair_data={}
        active=set()
        omitted=[]
        entered=0
        retained=0

        for i,j in filtered:
            edge=_uid_edge(
                uids[i],uids[j]
            )
            data=molecular_pair_data_v20(
                basis,i,j,self.provider,
                self.dt,qdots,pdots,
                centers,diag_h_abs,self.settings,
            )
            pair_data[(i,j)]=data
            threshold=(
                self.settings.exit_score
                if edge in old
                else self.settings.enter_score
            )
            if data.score>=threshold:
                active.add(edge)
                if edge in old:
                    retained+=1
                else:
                    entered+=1
            else:
                omitted.append(
                    (data.score,edge,(i,j))
                )

        # Global local omitted-score budget.
        omitted_sq=sum(
            float(score)**2
            for score,_,_ in omitted
        )
        promoted=0
        budget=float(
            self.settings.local_omitted_score_l2_budget
        )
        if math.sqrt(omitted_sq)>budget:
            for score,edge,ij in sorted(
                omitted,
                key=lambda x:(-x[0],x[1]),
            ):
                if math.sqrt(max(omitted_sq,0.0))<=budget:
                    break
                active.add(edge)
                omitted_sq-=float(score)**2
                promoted+=1
                if edge in old:
                    retained+=1
                else:
                    entered+=1

        self._active_uid_edges=active
        self.update_count+=1
        self.total_exact_pair_checks+=len(filtered)

        active_indices=tuple(sorted(
            (
                min(uid_to_index[a],uid_to_index[b]),
                max(uid_to_index[a],uid_to_index[b]),
            )
            for a,b in active
        ))

        active_index_set=set(active_indices)
        maximum_omitted=max(
            (
                data.score
                for ij,data in pair_data.items()
                if ij not in active_index_set
            ),
            default=0.0,
        )

        return MolecularSparseUpdateV20(
            active_edges=active_indices,
            pair_data=pair_data,
            diagonal_S=diag_S,
            diagonal_H=diag_H,
            diagonal_T=diag_T,
            qdots=qdots,
            pdots=pdots,
            spatial_candidate_pairs=spatial_count,
            candidate_pairs=len(filtered),
            exact_pair_checks=len(filtered),
            budget_promoted_edges=promoted,
            omitted_score_l2=float(
                math.sqrt(max(omitted_sq,0.0))
            ),
            total_offdiagonal_pairs=total,
            entered_edges=entered,
            exited_edges=len(old-active),
            retained_edges=retained,
            maximum_omitted_candidate_score=
                float(maximum_omitted),
        )


def build_sparse_molecular_matrices_v20(
    basis,
    update,
):
    n=len(basis)
    rowsS=[]; colsS=[]; dataS=[]
    rowsH=[]; colsH=[]; dataH=[]
    rowsT=[]; colsT=[]; dataT=[]

    for i in range(n):
        rowsS.append(i); colsS.append(i)
        dataS.append(update.diagonal_S[i])
        rowsH.append(i); colsH.append(i)
        dataH.append(update.diagonal_H[i])
        rowsT.append(i); colsT.append(i)
        dataT.append(update.diagonal_T[i])

    for i,j in update.active_edges:
        pair=update.pair_data[(i,j)]
        rowsS.extend((i,j)); colsS.extend((j,i))
        dataS.extend((pair.S,np.conj(pair.S)))

        rowsH.extend((i,j)); colsH.extend((j,i))
        dataH.extend((pair.H,np.conj(pair.H)))

        rowsT.extend((i,j)); colsT.extend((j,i))
        dataT.extend((pair.T_ij,pair.T_ji))

    shape=(n,n)
    return SparseMolecularMatricesV20(
        S=sparse.coo_matrix(
            (np.asarray(dataS,complex),(rowsS,colsS)),
            shape=shape,
        ).tocsr(),
        H=sparse.coo_matrix(
            (np.asarray(dataH,complex),(rowsH,colsH)),
            shape=shape,
        ).tocsr(),
        T_seed=sparse.coo_matrix(
            (np.asarray(dataT,complex),(rowsT,colsT)),
            shape=shape,
        ).tocsr(),
        active_edges=update.active_edges,
        update=update,
    )


def build_dense_molecular_reference_v20(
    basis,
    provider,
    dt,
    settings=MolecularSparseSettingsV20(),
):
    """Dense reference using the exact same v0.20 pair-centroid approximation."""
    basis=list(basis)
    qdots,pdots,centers=_center_kinematics(
        basis,provider
    )
    n=len(basis)

    diag_S=np.ones(n,dtype=complex)
    diag_H=np.zeros(n,dtype=complex)
    diag_T=np.zeros(n,dtype=complex)
    for i,b in enumerate(basis):
        point=centers[i].point
        diag_H[i]=(
            kinetic_matrix_element_equal_width(
                b.q,b.p,b.q,b.p,b.A,
                point.mass_matrix_q_au,
            )
            +point.energies[b.state]
        )
        diag_T[i]=(
            basis_time_matrix_element_equal_width(
                b.q,b.p,b.q,b.p,b.A,
                qdots[i],pdots[i],
            )
        )

    S=np.zeros((n,n),dtype=complex)
    H=np.zeros((n,n),dtype=complex)
    T=np.zeros((n,n),dtype=complex)
    np.fill_diagonal(S,diag_S)
    np.fill_diagonal(H,diag_H)
    np.fill_diagonal(T,diag_T)

    diag_h_abs=np.abs(diag_H)
    pair_data={}
    for i in range(n):
        for j in range(i+1,n):
            pair=molecular_pair_data_v20(
                basis,i,j,provider,float(dt),
                qdots,pdots,centers,
                diag_h_abs,settings,
            )
            pair_data[(i,j)]=pair
            S[i,j]=pair.S
            S[j,i]=np.conj(pair.S)
            H[i,j]=pair.H
            H[j,i]=np.conj(pair.H)
            T[i,j]=pair.T_ij
            T[j,i]=pair.T_ji

    return {
        "S":S,
        "H":H,
        "T_seed":T,
        "pair_data":pair_data,
        "pair_count":n*(n-1)//2,
    }


def relative_frobenius_error(A,B):
    A=np.asarray(A,complex)
    B=np.asarray(B,complex)
    return float(
        np.linalg.norm(A-B,ord="fro")
        /max(np.linalg.norm(B,ord="fro"),1e-30)
    )


def dense_audit_sparse_molecular_v20(
    basis,
    provider,
    dt,
    sparse_matrices,
    settings=MolecularSparseSettingsV20(),
):
    dense=build_dense_molecular_reference_v20(
        basis,provider,dt,settings
    )
    return {
        "relative_S_frobenius_error":
            relative_frobenius_error(
                sparse_matrices.S.toarray(),
                dense["S"],
            ),
        "relative_H_frobenius_error":
            relative_frobenius_error(
                sparse_matrices.H.toarray(),
                dense["H"],
            ),
        "relative_Tseed_frobenius_error":
            relative_frobenius_error(
                sparse_matrices.T_seed.toarray(),
                dense["T_seed"],
            ),
        "dense_pair_count":
            int(dense["pair_count"]),
        "active_edge_count":
            int(len(sparse_matrices.active_edges)),
    }
