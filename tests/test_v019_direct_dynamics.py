import numpy as np

from gaussian_dynamics.analytic_molecular_backend_v19 import (
    AnalyticMolecularLVCBackendV19,
    AnalyticMolecularLVCConfigV19,
    default_diatomic_two_mode_map_v19,
)
from gaussian_dynamics.molecular_direct_provider_v19 import (
    TrackedMolecularDirectProviderV19,
)
from gaussian_dynamics.local_gaussian_nd import (
    LocalAdiabaticTBF,
)
from gaussian_dynamics.molecular_direct_dynamics_v19 import (
    run_molecular_direct_dynamics_v19,
)


def _basis():
    A=0.8*np.eye(2)
    return [
        LocalAdiabaticTBF(
            0,np.array([-0.45,0.35]),
            np.array([8.0,0.0]),A
        ),
        LocalAdiabaticTBF(
            1,np.array([0.35,0.40]),
            np.array([-5.0,1.0]),A
        ),
    ]


def test_short_direct_dynamics_is_invariant_to_raw_root_scrambling():
    gmap=default_diatomic_two_mode_map_v19()
    clean=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
    )
    scrambled=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(
            gmap,
            AnalyticMolecularLVCConfigV19(
                scramble_roots=True
            ),
        ),
        gmap,
    )

    kwargs=dict(
        initial_basis=_basis(),
        C0=np.array([1.0+0j,0.2+0.1j]),
        dt=0.02,
        steps=20,
        spawn_threshold=1e9,
        max_basis=2,
        store_every=5,
    )
    a=run_molecular_direct_dynamics_v19(
        provider=clean,**kwargs
    )
    b=run_molecular_direct_dynamics_v19(
        provider=scrambled,**kwargs
    )

    assert np.allclose(
        a["norm"],b["norm"],atol=2e-10
    )
    assert np.allclose(
        a["final_coefficients"],
        b["final_coefficients"],
        atol=2e-9,
    )
    for x,y in zip(a["final_basis"],b["final_basis"]):
        assert np.allclose(x.q,y.q,atol=2e-10)
        assert np.allclose(x.p,y.p,atol=2e-10)

    assert not a["molecular_centroid_graph_audit"].get(
        "failed",False
    )
    assert (
        a["molecular_centroid_graph_audit"][
            "S_hermiticity_error"
        ]<1e-10
    )


def test_pyscf_v19_adapter_fails_cleanly_when_pyscf_missing():
    from gaussian_dynamics.pyscf_backend_v05 import (
        PySCFSACASSCFConfig,
    )
    from gaussian_dynamics.pyscf_molecular_bridge_v19 import (
        PySCFRawSnapshotBackendV19,
    )

    config=PySCFSACASSCFConfig(
        basis="sto-3g",
        ncas=2,
        nelecas=(1,1),
        nstates=2,
    )
    backend=PySCFRawSnapshotBackendV19(config)

    # Construction itself must not require PySCF.
    assert backend.engine.config.nstates==2


def test_overlap_engine_path_works_without_finite_state_vectors():
    from dataclasses import dataclass
    from gaussian_dynamics.molecular_backend import (
        CartesianElectronicStructurePoint,
    )
    from gaussian_dynamics.molecular_snapshot_v19 import (
        MolecularElectronicSnapshotV19,
    )

    @dataclass
    class FakeWave:
        vectors: np.ndarray

        def with_transformed_roots(self,permutation,phase_factors):
            perm=np.asarray(permutation,dtype=int)
            phase=np.asarray(phase_factors,dtype=complex)
            return FakeWave(
                self.vectors[:,perm]*phase[None,:]
            )

    class FakeBackend:
        def __init__(self,gmap):
            self.gmap=gmap
            self.clean=AnalyticMolecularLVCBackendV19(
                gmap,
                AnalyticMolecularLVCConfigV19(
                    scramble_roots=True
                ),
            )

        def evaluate_snapshot(self,geometry):
            raw=self.clean.evaluate_snapshot(geometry)
            return MolecularElectronicSnapshotV19(
                point=raw.point,
                wavefunction_snapshot=FakeWave(
                    raw.state_vectors
                ),
            ).validate()

    def overlap_engine(previous,current):
        return (
            previous.wavefunction_snapshot.vectors.conj().T
            @current.wavefunction_snapshot.vectors
        )

    gmap=default_diatomic_two_mode_map_v19()
    provider=TrackedMolecularDirectProviderV19(
        FakeBackend(gmap),
        gmap,
        overlap_engine=overlap_engine,
    )
    clean=TrackedMolecularDirectProviderV19(
        AnalyticMolecularLVCBackendV19(gmap),
        gmap,
    )

    for q in (
        np.array([-0.5,0.35]),
        np.array([-0.1,0.35]),
        np.array([0.3,0.35]),
    ):
        a=provider.evaluate(q)
        b=clean.evaluate(q)
        assert np.allclose(a.energies,b.energies,atol=1e-12)
        assert np.allclose(a.nac_q,b.nac_q,atol=1e-12)
