import numpy as np
import pytest

from gaussian_dynamics.adiabatic import adiabatic_point
from gaussian_dynamics.electronic_structure import (
    AnalyticAvoidedCrossingProvider,
    TabulatedElectronicStructureProvider,
    CachedProvider,
    project_cartesian_vector_to_coordinate,
    point_fingerprint,
)
from gaussian_dynamics.provider_dynamics import (
    ProviderTBF, velocity_verlet_tbf, energy_conserving_child
)


def test_analytic_provider_matches_v02():
    provider=AnalyticAvoidedCrossingProvider()
    q=0.37
    p=provider.evaluate(q)
    E,U,g,d=adiabatic_point(q)
    assert np.allclose(p.energies,E)
    assert np.allclose(p.gradients_q,g)
    assert np.allclose(p.nac_q,d)


def test_tabulated_provider_reproduces_nodes_and_interpolates():
    base=AnalyticAvoidedCrossingProvider()
    q=np.linspace(-1,1,9)
    pts=[base.evaluate(x) for x in q]
    tab=TabulatedElectronicStructureProvider(
        q,
        np.array([p.energies for p in pts]),
        np.array([p.gradients_q for p in pts]),
        np.array([p.nac_q for p in pts]),
    )
    for qi,ref in zip(q,pts):
        got=tab.evaluate(qi)
        assert np.allclose(got.energies,ref.energies)
        assert np.allclose(got.gradients_q,ref.gradients_q)
        assert np.allclose(got.nac_q,ref.nac_q)


def test_projection_chain_rule():
    grad=np.array([[1.,2.,3.],[-1.,0.5,4.]])
    tangent=np.array([[0.1,0.0,0.2],[0.0,-0.4,0.3]])
    expected=np.sum(grad*tangent)
    assert project_cartesian_vector_to_coordinate(grad,tangent) == pytest.approx(expected)


def test_cache_and_fingerprint_are_deterministic():
    wrapped=CachedProvider(AnalyticAvoidedCrossingProvider())
    a=wrapped.evaluate(0.123456789)
    b=wrapped.evaluate(0.123456789)
    assert wrapped.misses == 1
    assert wrapped.hits == 1
    assert point_fingerprint(a) == point_fingerprint(b)


def test_provider_velocity_verlet_is_finite():
    provider=CachedProvider(AnalyticAvoidedCrossingProvider())
    tbf=ProviderTBF(state=0,q=-1.0,p=0.4,alpha=1.0)
    new=velocity_verlet_tbf(tbf,provider,mass=20.0,dt=0.01)
    assert np.isfinite(new.q)
    assert np.isfinite(new.p)


def test_provider_child_conserves_energy_when_allowed():
    provider=AnalyticAvoidedCrossingProvider()
    parent=ProviderTBF(state=1,q=0.0,p=1.0)
    child=energy_conserving_child(parent,provider,mass=20.0)
    assert child is not None
    p=provider.evaluate(parent.q)
    before=parent.p**2/(40.0)+p.energies[parent.state]
    after=child.p**2/(40.0)+p.energies[child.state]
    assert abs(before-after) < 1e-12


def test_pyscf_provider_has_clear_optional_dependency_behavior():
    from gaussian_dynamics.pyscf_provider import PySCFStateAveragedCASSCFProvider

    def builder(q):
        return [("H",(0,0,-q/2)),("H",(0,0,q/2))], np.array([[0,0,-.5],[0,0,.5]])

    provider=PySCFStateAveragedCASSCFProvider(
        builder,"sto-3g",2,2,nstates=2,verbose=0
    )

    try:
        import pyscf  # noqa
    except ImportError:
        with pytest.raises(ImportError, match="PySCF is optional"):
            provider.evaluate(1.5)
