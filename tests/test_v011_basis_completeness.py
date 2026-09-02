import numpy as np

from gaussian_dynamics.basis_completeness import (
    overlap_spectrum_metrics,
    canonical_coefficient_weights,
    width_diversity_ratio,
    generation_histogram,
)


class Dummy:
    def __init__(self,uid,A):
        self.uid=uid
        self.A=np.asarray(A,float)


def test_overlap_spectrum_and_canonical_weights():
    S=np.array([
        [1.0,0.2],
        [0.2,1.0],
    ],complex)
    C=np.array([1.0,0.3j])

    report=overlap_spectrum_metrics(S)
    weights=canonical_coefficient_weights(C,S)

    norm=float(np.real(np.vdot(C,S@C)))

    assert report["numerical_rank"]==2
    assert report["condition_number"]>1.0
    assert np.isclose(weights["norm"],norm)
    assert 1.0 <= weights["participation_ratio"] <= 2.0


def test_width_diversity_and_generation_histogram():
    basis=[
        Dummy(0,np.eye(2)),
        Dummy(1,2*np.eye(2)),
        Dummy(2,0.5*np.eye(2)),
    ]
    lineage={
        0:{"generation":0},
        1:{"generation":1},
        2:{"generation":1},
    }

    assert np.isclose(width_diversity_ratio(basis),16.0)
    assert generation_histogram(lineage,basis)=={0:1,1:2}
