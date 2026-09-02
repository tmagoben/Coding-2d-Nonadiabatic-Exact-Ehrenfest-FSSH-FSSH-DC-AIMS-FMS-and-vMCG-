from dataclasses import dataclass
import numpy as np

from .ci2d import (
    LVC2DParameters,
    adiabatic_energies_2d,
    adiabatic_gradients_2d,
    vector_nac_2d,
    analytic_adiabatic_vectors,
)
from .molecular_backend import (
    AMU_TO_ELECTRON_MASS,
    LinearGeometryMap,
    MolecularGeometry,
    CartesianElectronicStructurePoint,
)
from .molecular_snapshot_v19 import MolecularElectronicSnapshotV19


def default_diatomic_two_mode_map_v19(
    bond_length_bohr=1.4,
):
    """Synthetic H2 geometry with two orthonormal Cartesian collective modes."""
    R0=np.array([
        [-0.5*bond_length_bohr,0.0,0.0],
        [ 0.5*bond_length_bohr,0.0,0.0],
    ])
    inv=np.sqrt(0.5)
    modes=np.zeros((2,2,3),dtype=float)

    # Stretch-like relative x displacement.
    modes[0,0,0]=-inv
    modes[0,1,0]= inv

    # Relative transverse y displacement.
    modes[1,0,1]= inv
    modes[1,1,1]=-inv

    return LinearGeometryMap(("H","H"),R0,modes)


@dataclass(frozen=True)
class AnalyticMolecularLVCConfigV19:
    params: LVC2DParameters=LVC2DParameters()
    masses_amu: tuple=(1.00784,1.00784)
    scramble_roots: bool=False
    fail_if_q0_greater_than: float | None=None


class AnalyticMolecularLVCBackendV19:
    """Deterministic molecular-style backend with Cartesian gradients/NACs.

    It embeds the exact 2D LVC electronic problem into a two-atom Cartesian geometry.
    The backend may deliberately swap/sign-flip raw roots to validate state tracking.
    """

    def __init__(
        self,
        geometry_map=None,
        config=AnalyticMolecularLVCConfigV19(),
    ):
        self.geometry_map=(
            geometry_map
            if geometry_map is not None
            else default_diatomic_two_mode_map_v19()
        )
        if self.geometry_map.nq!=2:
            raise ValueError("AnalyticMolecularLVCBackendV19 requires exactly two modes.")
        self.config=config
        self.calls=0

        J=self.geometry_map.J
        self._pinv=np.linalg.pinv(J)
        self._cart_lift=J@np.linalg.inv(J.T@J)

    def generalized_coordinates(self,geometry):
        geometry=MolecularGeometry(
            geometry.symbols,geometry.coords_bohr
        )
        if geometry.symbols!=self.geometry_map.symbols:
            raise ValueError("geometry symbols differ from the configured map.")
        dr=(geometry.coords_bohr-self.geometry_map.reference_bohr).reshape(-1)
        return self._pinv@dr

    def _raw_transform(self,q):
        perm=np.array([0,1],dtype=int)
        phase=np.ones(2,dtype=float)

        if self.config.scramble_roots:
            # Deliberately discontinuous raw root order and signs.
            if q[0]>0.0:
                perm=np.array([1,0],dtype=int)
            if q[1]<0.0:
                phase[0]*=-1.0
            if q[0]+0.5*q[1]>0.35:
                phase[1]*=-1.0
        return perm,phase

    def evaluate_snapshot(self,geometry):
        self.calls+=1
        q=np.asarray(
            self.generalized_coordinates(geometry),
            dtype=float,
        )

        if (
            self.config.fail_if_q0_greater_than is not None
            and q[0]>float(self.config.fail_if_q0_greater_than)
        ):
            raise RuntimeError(
                "Deterministic v0.19 backend failure region reached."
            )

        E=adiabatic_energies_2d(
            q[0],q[1],self.config.params
        )
        Gq=adiabatic_gradients_2d(
            q,self.config.params
        )
        Dq=vector_nac_2d(
            q,self.config.params
        )
        V=analytic_adiabatic_vectors(
            q,self.config.params
        )

        # Lift generalized covectors into Cartesian space so J^T g_R = g_q.
        Gcart=(self._cart_lift@Gq.T).T.reshape(2,2,3)
        Dflat=np.einsum(
            "ra,ija->ijr",
            self._cart_lift,Dq,
        )
        Dcart=Dflat.reshape(2,2,2,3)

        perm,phase=self._raw_transform(q)

        Eraw=E[perm]
        Graw=Gcart[perm]
        Draw=Dcart[np.ix_(perm,perm)].copy()
        phase_matrix=phase[:,None]*phase[None,:]
        Draw=Draw*phase_matrix[:,:,None,None]
        Vraw=V[:,perm]*phase[None,:]

        masses=np.asarray(
            self.config.masses_amu,
            dtype=float,
        )
        point=CartesianElectronicStructurePoint(
            geometry=geometry,
            energies=Eraw,
            gradients_cart=Graw,
            nac_cart=Draw,
            masses_amu=masses,
            metadata={
                "backend":"analytic_molecular_lvc_v19",
                "raw_permutation":perm.tolist(),
                "raw_phase":phase.tolist(),
                "q_reconstructed":q.tolist(),
                "scramble_roots":bool(self.config.scramble_roots),
                "dynamics_mult_ediff":False,
            },
        ).validate()

        return MolecularElectronicSnapshotV19(
            point=point,
            state_vectors=Vraw,
            metadata={
                "exact_generalized_coordinates":q.tolist(),
            },
        ).validate()

    def evaluate(self,geometry):
        return self.evaluate_snapshot(geometry).point

    @property
    def generalized_mass_matrix_au(self):
        return self.geometry_map.mass_matrix_q_au(
            np.asarray(self.config.masses_amu,float)
        )
