from dataclasses import dataclass, field
import numpy as np

from .molecular_backend import (
    CartesianElectronicStructurePoint,
    GeneralizedElectronicStructurePoint,
)


@dataclass
class MolecularElectronicSnapshotV19:
    """Electronic point plus a cross-geometry state representation.

    `state_vectors[:,i]` is the raw electronic ket for state i in a common validation
    representation.  Real PySCF/CASSCF use should instead supply a wavefunction
    snapshot and an injected overlap engine; the finite vectors are primarily for the
    deterministic v0.19 molecular benchmark.
    """
    point: CartesianElectronicStructurePoint
    state_vectors: np.ndarray | None = None
    wavefunction_snapshot: object | None = None
    metadata: dict = field(default_factory=dict)

    def validate(self):
        self.point=self.point.validate()
        ns=len(self.point.energies)

        if self.state_vectors is not None:
            V=np.asarray(self.state_vectors,dtype=complex)
            if V.ndim!=2 or V.shape[1]!=ns:
                raise ValueError(
                    "state_vectors must have shape (representation_dimension,nstate)."
                )
            gram=V.conj().T@V
            if not np.allclose(gram,np.eye(ns),atol=1e-10):
                raise ValueError("state_vectors must be orthonormal.")
            self.state_vectors=V

        if self.state_vectors is None and self.wavefunction_snapshot is None:
            raise ValueError(
                "A molecular snapshot requires state_vectors or wavefunction_snapshot."
            )
        return self


@dataclass
class TrackedGeneralizedSnapshotV19:
    point: GeneralizedElectronicStructurePoint
    state_vectors: np.ndarray | None
    wavefunction_snapshot: object | None
    node_id: str
    tracking_metadata: dict
    source_metadata: dict = field(default_factory=dict)

    def validate(self):
        self.point=self.point.validate()
        if self.state_vectors is not None:
            V=np.asarray(self.state_vectors,dtype=complex)
            ns=len(self.point.energies)
            if V.ndim!=2 or V.shape[1]!=ns:
                raise ValueError("tracked state_vectors have incompatible shape.")
            if not np.allclose(V.conj().T@V,np.eye(ns),atol=1e-10):
                raise ValueError("tracked state_vectors must be orthonormal.")
            self.state_vectors=V
        return self
