import numpy as np
import pytest

from gaussian_dynamics.complex_operator_cache_v213 import (
    FixedFrameComplexOperatorCacheV213,
)
from gaussian_dynamics.electronic_contract_v213 import (
    ElectronicModelSpaceV213,
    ElectronicOperatorProvenanceV213,
    ElectronicStateDescriptorV213,
)
from gaussian_dynamics.synthetic_operator_provider_v21 import (
    SyntheticLinearOperatorConfigV21,
    SyntheticLinearOperatorProviderV21,
)
from gaussian_dynamics.electronic_operator_v21 import ElectronicOperatorSnapshotV21


def _provenance(parameter):
    space = ElectronicModelSpaceV213(
        name="three-state fixed fixture",
        representation="fixed_general",
        states=tuple(
            ElectronicStateDescriptorV213(f"state-{index}") for index in range(3)
        ),
    )
    return ElectronicOperatorProvenanceV213(
        model_name="synthetic-linear",
        model_version="1",
        model_space=space,
        spin_free_method="analytic fixture",
        parameters={"parameter": parameter},
    )


def test_complex_cache_roundtrip_and_provenance_separation(tmp_path):
    base = SyntheticLinearOperatorProviderV21(
        SyntheticLinearOperatorConfigV21(nstate=3, nq=2, seed=21301)
    )
    first = FixedFrameComplexOperatorCacheV213(
        base, tmp_path, _provenance(1.0), namespace="contract-test"
    )
    q = np.asarray([0.17, -0.28])
    miss = first.evaluate_snapshot(q)
    hit = first.evaluate_snapshot(q)

    assert base.calls == 1
    assert first.misses == 1
    assert first.hits == 1
    assert np.array_equal(hit.point.H, miss.point.H)
    assert np.array_equal(hit.point.dH_dq, miss.point.dH_dq)
    assert np.iscomplexobj(hit.point.H)
    assert hit.metadata["cache_status"] == "hit"

    second = FixedFrameComplexOperatorCacheV213(
        base, tmp_path, _provenance(2.0), namespace="contract-test"
    )
    other = second.evaluate_snapshot(q)
    assert base.calls == 2
    assert second.misses == 1
    assert first.provider_fingerprint != second.provider_fingerprint
    assert np.array_equal(other.point.H, miss.point.H)
    assert len(tuple(tmp_path.glob("*.npz"))) == 2
    assert len(tuple(tmp_path.glob("*.json"))) == 2

    rotated = ElectronicOperatorSnapshotV21(
        point=hit.point,
        state_vectors=np.diag([1.0, 1.0j, -1.0]).astype(complex),
    ).validate()
    with pytest.raises(ValueError, match="fixed-frame cross-geometry overlap"):
        first.snapshot_overlap(hit, rotated)

    class _WavefunctionBearingProvider:
        def evaluate_snapshot(self, coordinates):
            snapshot = base.evaluate_snapshot(coordinates)
            snapshot.wavefunction_snapshot = object()
            return snapshot

    refusing = FixedFrameComplexOperatorCacheV213(
        _WavefunctionBearingProvider(),
        tmp_path / "moving-refused",
        _provenance(3.0),
    )
    with pytest.raises(ValueError, match="dedicated molecular cache"):
        refusing.evaluate_snapshot(q)
