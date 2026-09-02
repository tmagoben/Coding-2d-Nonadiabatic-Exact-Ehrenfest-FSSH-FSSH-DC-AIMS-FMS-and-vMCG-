import numpy as np

from gaussian_dynamics.state_tracking import (
    maximum_overlap_assignment,
)
from gaussian_dynamics.state_tracking_v19 import (
    scalable_maximum_overlap_assignment_v19,
)


def _random_unitary(n,seed):
    rng=np.random.default_rng(seed)
    X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,R=np.linalg.qr(X)
    phase=np.diag(R)
    phase=np.where(np.abs(phase)>0,phase/np.abs(phase),1.0)
    return Q*phase.conj()[None,:]


def test_hungarian_tracking_matches_exhaustive_small_manifolds():
    for n in (2,3,4,5):
        O=_random_unitary(n,100+n)
        a=maximum_overlap_assignment(
            O,
            minimum_overlap=0.0,
            minimum_score_margin=0.0,
            real_gauge=False,
        )
        b=scalable_maximum_overlap_assignment_v19(
            O,
            minimum_overlap=0.0,
            minimum_score_margin=0.0,
            real_gauge=False,
        )
        assert np.array_equal(a.permutation,b.permutation)
        assert np.isclose(a.best_score,b.best_score,atol=1e-12)
        assert np.isclose(
            a.second_best_score,b.second_best_score,atol=1e-12
        )
        assert np.isclose(a.score_margin,b.score_margin,atol=1e-12)


def test_hungarian_tracking_handles_eight_states_without_factorial_search():
    O=_random_unitary(8,777)
    out=scalable_maximum_overlap_assignment_v19(
        O,
        minimum_overlap=0.0,
        minimum_score_margin=0.0,
        real_gauge=False,
    )
    assert sorted(out.permutation.tolist())==list(range(8))
    assert np.isfinite(out.best_score)
    assert np.isfinite(out.second_best_score)
