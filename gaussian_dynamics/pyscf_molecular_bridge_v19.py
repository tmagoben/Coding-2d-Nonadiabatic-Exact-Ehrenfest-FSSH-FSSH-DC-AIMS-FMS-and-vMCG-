from .molecular_snapshot_v19 import MolecularElectronicSnapshotV19
from .pyscf_tracked_backend_v06 import (
    PySCFTrackedSACASSCFBackend,
)
from .pyscf_wavefunction_overlap import (
    casscf_state_overlap_matrix,
)


class PySCFRawSnapshotBackendV19:
    """Raw SA-CASSCF point/snapshot adapter for the v0.19 molecular provider.

    State tracking is deliberately *not* performed here. The nearest-anchor v0.19
    provider owns state identity for branched TBF centers/pair centroids. This avoids
    call-order-dependent sequential tracking.

    PySCF is imported only when a calculation is actually evaluated.
    """

    def __init__(self,config):
        self.engine=PySCFTrackedSACASSCFBackend(
            config,
            ambiguity_policy="accept",
        )

    def evaluate_snapshot(self,geometry):
        point,snapshot=(
            self.engine.evaluate_raw_with_snapshot(
                geometry
            )
        )
        return MolecularElectronicSnapshotV19(
            point=point,
            wavefunction_snapshot=snapshot,
            metadata={
                "backend_adapter":
                    "PySCFRawSnapshotBackendV19",
                "state_tracking":
                    "deferred_to_v19_nearest_anchor_provider",
            },
        ).validate()

    def evaluate(self,geometry):
        return self.evaluate_snapshot(geometry).point


def pyscf_snapshot_overlap_engine_v19(
    previous,
    current,
):
    """Cross-geometry many-electron overlap for v0.19 tracked snapshots."""
    return casscf_state_overlap_matrix(
        previous.wavefunction_snapshot,
        current.wavefunction_snapshot,
    )
