from dataclasses import dataclass
import numpy as np

from .gaussian_general import validate_spd
from .lvc_exact_gaussian import IDENTITY_2, SIGMA_X, SIGMA_Z


@dataclass(frozen=True)
class GaussianPairData:
    """Reusable algebra for one oriented cross Gaussian g_i^* g_j.

    A single dense solve of

        B = A_i + A_j

    supplies both

        mu = B^{-1} l
        Sigma = B^{-1}.

    The same overlap, cross centroid, and covariance are then reused by overlap,
    kinetic, LVC potential, and moving-basis time matrix elements.
    """

    qi: np.ndarray
    pi: np.ndarray
    Ai: np.ndarray
    qj: np.ndarray
    pj: np.ndarray
    Aj: np.ndarray

    overlap: complex
    centroid: np.ndarray
    covariance: np.ndarray

    def reversed(self):
        """Return the j,i orientation without another dense pair solve."""
        return GaussianPairData(
            qi=self.qj,
            pi=self.pj,
            Ai=self.Aj,
            qj=self.qi,
            pj=self.pi,
            Aj=self.Ai,
            overlap=np.conj(self.overlap),
            centroid=np.conj(self.centroid),
            covariance=self.covariance,
        )

    def kinetic(self, mass_inverse):
        """Exact nuclear kinetic matrix element using cached pair moments."""
        Minv=np.asarray(mass_inverse,dtype=float)

        ui=-self.Ai@(self.centroid-self.qi)-1j*self.pi
        uj=-self.Aj@(self.centroid-self.qj)+1j*self.pj

        fluctuation=np.trace(
            self.Ai@Minv@self.Aj@self.covariance
        )
        return 0.5*self.overlap*(
            ui@Minv@uj+fluctuation
        )

    def lvc_potential_matrix(self, params):
        """Exact 2x2 LVC potential matrix element from cached moments."""
        mu=self.centroid
        Sigma=self.covariance

        second_sum=(
            mu[0]*mu[0]
            +mu[1]*mu[1]
            +Sigma[0,0]
            +Sigma[1,1]
        )
        common=0.5*params.omega**2*second_sum

        return self.overlap*(
            common*IDENTITY_2
            +params.kappa*mu[0]*SIGMA_Z
            +params.lam*mu[1]*SIGMA_X
        )

    def time_element(
        self,
        qdot_j,
        pdot_j,
        Adot_j=None,
    ):
        """Exact <g_i|dot g_j> reusing the same pair moments."""
        qdot_j=np.asarray(qdot_j,dtype=float)
        pdot_j=np.asarray(pdot_j,dtype=float)

        if Adot_j is None:
            Adot_j=np.zeros_like(self.Aj)
        Adot_j=np.asarray(Adot_j,dtype=float)
        if (
            Adot_j.shape!=self.Aj.shape
            or not np.allclose(Adot_j,Adot_j.T,atol=1e-12)
        ):
            raise ValueError(
                "Adot_j must be symmetric and match Aj."
            )

        y=self.centroid-self.qj

        normalization=0.25*np.trace(
            np.linalg.solve(self.Aj,Adot_j)
        )
        center=(self.Aj@y-1j*self.pj)@qdot_j
        momentum=1j*y@pdot_j
        width=-0.5*(
            y@Adot_j@y
            +np.trace(Adot_j@self.covariance)
        )

        return self.overlap*(
            normalization+center+momentum+width
        )


@dataclass
class PairCacheStats:
    requests: int = 0
    canonical_solves: int = 0
    direct_hits: int = 0
    reverse_views: int = 0
    inherited_pairs: int = 0

    def as_dict(self):
        return {
            "requests":int(self.requests),
            "canonical_solves":int(self.canonical_solves),
            "direct_hits":int(self.direct_hits),
            "reverse_views":int(self.reverse_views),
            "inherited_pairs":int(self.inherited_pairs),
        }


def _state_arrays(tbf):
    q=np.asarray(tbf.q,dtype=float)
    p=np.asarray(tbf.p,dtype=float)
    A=validate_spd(tbf.A)

    if q.ndim!=1 or p.shape!=q.shape:
        raise ValueError("TBF q and p must be equal-length vectors.")
    if A.shape!=(len(q),len(q)):
        raise ValueError("TBF width has incompatible shape.")
    return q,p,A


def _build_pair_data(
    tbf_i,
    tbf_j,
    logdet_Ai=None,
    logdet_Aj=None,
):
    qi,pi,Ai=_state_arrays(tbf_i)
    qj,pj,Aj=_state_arrays(tbf_j)

    if qi.shape!=qj.shape:
        raise ValueError("pair TBF dimensions do not match.")

    D=len(qi)
    B=Ai+Aj
    ell=Ai@qi+Aj@qj+1j*(pj-pi)

    # One dense solve/factorization for BOTH mu and Sigma.
    rhs=np.column_stack([
        ell,
        np.eye(D,dtype=complex),
    ])
    solved=np.linalg.solve(
        B.astype(complex),
        rhs,
    )
    mu=solved[:,0]
    Sigma=np.real_if_close(
        solved[:,1:],
        tol=1000,
    ).astype(float)

    if logdet_Ai is None:
        logdet_Ai=np.linalg.slogdet(Ai)[1]
    if logdet_Aj is None:
        logdet_Aj=np.linalg.slogdet(Aj)[1]
    logdet_B=np.linalg.slogdet(B)[1]

    c=(
        -0.5*qi@Ai@qi
        -0.5*qj@Aj@qj
        +1j*pi@qi
        -1j*pj@qj
    )
    log_prefactor=(
        0.25*logdet_Ai
        +0.25*logdet_Aj
        -0.5*logdet_B
        +0.5*D*np.log(2.0)
    )
    overlap=np.exp(
        log_prefactor+c+0.5*ell@mu
    )

    return GaussianPairData(
        qi=qi,
        pi=pi,
        Ai=Ai,
        qj=qj,
        pj=pj,
        Aj=Aj,
        overlap=complex(overlap),
        centroid=np.asarray(mu,dtype=complex),
        covariance=np.asarray(Sigma,dtype=float),
    )


class GaussianPairCache:
    """Lazy canonical-pair cache for one frozen Gaussian basis snapshot.

    Canonical pairs use i <= j.  A reversed orientation is obtained by exact
    conjugation/swap identities and therefore requires no additional B=A_i+A_j solve.

    The cache is valid only while the contained basis snapshot is unchanged.
    """

    def __init__(
        self,
        basis,
        *,
        inherited=None,
        inherited_index_map=None,
    ):
        self.basis=list(basis)
        if not self.basis:
            raise ValueError("pair cache basis cannot be empty.")

        self.dimension=len(np.asarray(self.basis[0].q,float))
        self._logdet_A=[
            float(np.linalg.slogdet(validate_spd(b.A))[1])
            for b in self.basis
        ]
        self._pairs={}
        self.stats=PairCacheStats()

        if inherited is not None:
            if inherited_index_map is None:
                raise ValueError(
                    "inherited_index_map is required with inherited cache."
                )
            index_map={
                int(new):int(old)
                for new,old in inherited_index_map.items()
            }
            old_to_new={
                old:new
                for new,old in index_map.items()
            }

            for (oi,oj),data in inherited._pairs.items():
                if oi not in old_to_new or oj not in old_to_new:
                    continue
                ni=old_to_new[oi]
                nj=old_to_new[oj]
                if ni<=nj:
                    self._pairs[(ni,nj)]=data
                else:
                    self._pairs[(nj,ni)]=data.reversed()
                self.stats.inherited_pairs+=1

    def __len__(self):
        return len(self.basis)

    def pair(self,i,j):
        i=int(i)
        j=int(j)
        n=len(self.basis)
        if not (0<=i<n and 0<=j<n):
            raise IndexError("pair index out of range.")

        self.stats.requests+=1

        if i<=j:
            key=(i,j)
            reverse=False
        else:
            key=(j,i)
            reverse=True

        if key in self._pairs:
            self.stats.direct_hits+=1
            data=self._pairs[key]
        else:
            a,b=key
            data=_build_pair_data(
                self.basis[a],
                self.basis[b],
                self._logdet_A[a],
                self._logdet_A[b],
            )
            self._pairs[key]=data
            self.stats.canonical_solves+=1

        if reverse:
            self.stats.reverse_views+=1
            return data.reversed()
        return data

    def prime(self):
        """Build every canonical pair exactly once."""
        for i in range(len(self.basis)):
            for j in range(i,len(self.basis)):
                self.pair(i,j)
        return self

    @property
    def canonical_pair_count(self):
        n=len(self.basis)
        return n*(n+1)//2

    @property
    def cached_pair_count(self):
        return len(self._pairs)

    def expanded(self,child):
        """New snapshot cache reusing all old-old pair data at the same geometry."""
        new_basis=self.basis+[child]
        mapping={i:i for i in range(len(self.basis))}
        return GaussianPairCache(
            new_basis,
            inherited=self,
            inherited_index_map=mapping,
        )

    def subset(self,keep):
        """Subset cache without recomputing surviving pair data."""
        keep=[int(i) for i in keep]
        new_basis=[self.basis[i] for i in keep]
        mapping={
            new:old
            for new,old in enumerate(keep)
        }
        return GaussianPairCache(
            new_basis,
            inherited=self,
            inherited_index_map=mapping,
        )


def build_cached_spinor_lvc_matrices(
    cache,
    provider,
):
    """Build exact spinor-complete S/H using one cached pair solve per i<=j."""
    n=len(cache)
    ns=2
    dim=ns*n

    S=np.zeros((dim,dim),dtype=complex)
    H=np.zeros((dim,dim),dtype=complex)
    Snuc=np.zeros((n,n),dtype=complex)

    point=provider.evaluate(
        np.asarray(cache.basis[0].q,float)
    )
    Minv=np.linalg.inv(
        np.asarray(point.mass_matrix,float)
    )
    params=provider.params
    eye=np.eye(ns,dtype=complex)

    for i in range(n):
        si=slice(ns*i,ns*(i+1))
        for j in range(i,n):
            sj=slice(ns*j,ns*(j+1))
            pair=cache.pair(i,j)

            Sij=pair.overlap
            Tij=pair.kinetic(Minv)
            Vij=pair.lvc_potential_matrix(params)

            blockS=Sij*eye
            blockH=Tij*eye+Vij

            Snuc[i,j]=Sij
            S[si,sj]=blockS
            H[si,sj]=blockH

            if i!=j:
                Snuc[j,i]=np.conj(Sij)
                S[sj,si]=blockS.conj().T
                H[sj,si]=blockH.conj().T

    return S,H,Snuc


def build_cached_spinor_time_matrix(
    cache,
    qdots,
    pdots,
    Adots=None,
):
    """Build ordered T_ij while reusing canonical pair algebra.

    T is not Hermitian, so all ordered i,j entries are evaluated.  The expensive
    cross-Gaussian B solve is nevertheless performed only once per canonical pair.
    """
    n=len(cache)
    ns=2
    dim=ns*n

    qdots=np.asarray(qdots,float)
    pdots=np.asarray(pdots,float)
    nq=cache.dimension

    if qdots.shape!=(n,nq) or pdots.shape!=qdots.shape:
        raise ValueError("qdots/pdots have incompatible shape.")

    if Adots is None:
        Adots=[None]*n
    if len(Adots)!=n:
        raise ValueError("Adots must contain one entry per Gaussian.")

    T=np.zeros((dim,dim),dtype=complex)
    eye=np.eye(ns,dtype=complex)

    for i in range(n):
        si=slice(ns*i,ns*(i+1))
        for j in range(n):
            sj=slice(ns*j,ns*(j+1))
            pair=cache.pair(i,j)
            tij=pair.time_element(
                qdots[j],
                pdots[j],
                Adots[j],
            )
            T[si,sj]=tij*eye

    return T


def expand_cached_spinor_lvc_matrices(
    S_old,
    H_old,
    Snuc_old,
    expanded_cache,
    provider,
):
    """Add the final Gaussian row/column without rebuilding old-old blocks."""
    nnew=len(expanded_cache)
    nold=nnew-1
    ns=2

    S0=np.asarray(S_old,dtype=complex)
    H0=np.asarray(H_old,dtype=complex)
    N0=np.asarray(Snuc_old,dtype=complex)

    if S0.shape!=(ns*nold,ns*nold):
        raise ValueError("S_old has incompatible dimension.")
    if H0.shape!=S0.shape or N0.shape!=(nold,nold):
        raise ValueError("old matrices have incompatible shapes.")

    S=np.zeros((ns*nnew,ns*nnew),dtype=complex)
    H=np.zeros_like(S)
    Snuc=np.zeros((nnew,nnew),dtype=complex)

    S[:ns*nold,:ns*nold]=S0
    H[:ns*nold,:ns*nold]=H0
    Snuc[:nold,:nold]=N0

    point=provider.evaluate(
        np.asarray(expanded_cache.basis[0].q,float)
    )
    Minv=np.linalg.inv(
        np.asarray(point.mass_matrix,float)
    )
    params=provider.params
    eye=np.eye(ns,dtype=complex)

    j=nold
    sj=slice(ns*j,ns*(j+1))

    for i in range(nnew):
        si=slice(ns*i,ns*(i+1))
        pair=expanded_cache.pair(i,j)

        Sij=pair.overlap
        Tij=pair.kinetic(Minv)
        Vij=pair.lvc_potential_matrix(params)

        blockS=Sij*eye
        blockH=Tij*eye+Vij

        Snuc[i,j]=Sij
        S[si,sj]=blockS
        H[si,sj]=blockH

        if i!=j:
            Snuc[j,i]=np.conj(Sij)
            S[sj,si]=blockS.conj().T
            H[sj,si]=blockH.conj().T

    return S,H,Snuc


def subset_cached_spinor_lvc_matrices(
    S,
    H,
    Snuc,
    cache,
    keep,
):
    """Remove Gaussian blocks by exact slicing; no pair algebra is recomputed."""
    keep=np.asarray(keep,dtype=int)
    ns=2

    eidx=np.array([
        ns*i+a
        for i in keep
        for a in range(ns)
    ],dtype=int)

    return (
        np.asarray(S)[np.ix_(eidx,eidx)].copy(),
        np.asarray(H)[np.ix_(eidx,eidx)].copy(),
        np.asarray(Snuc)[np.ix_(keep,keep)].copy(),
        cache.subset(keep),
    )


def v14_factorization_equivalent_for_sh(n_basis):
    """Dense B-solve/inverse equivalents used by the v0.14 S/H implementation.

    Per canonical pair:
      overlap helper                         1
      kinetic: overlap + centroid + Sigma   3
      LVC V: overlap + centroid + Sigma     3
    total                                   7
    """
    n=int(n_basis)
    return 7*n*(n+1)//2


def v14_factorization_equivalent_for_time(n_basis):
    """v0.14 moving-basis T: 3 pair solves/inverses for every ordered pair."""
    n=int(n_basis)
    return 3*n*n


def v15_factorization_count_for_snapshot(n_basis):
    """v0.15 cache: one multi-RHS solve per canonical pair snapshot."""
    n=int(n_basis)
    return n*(n+1)//2
