"""Independent NAC integration, finite differences, and gauge-covariant tangents."""
from dataclasses import replace

import numpy as np
import pytest

from gaussian_dynamics.analytic_soc_bridge_v307 import AnalyticSOCFrameBridgeV307, analytic_soc_fixture_v307
from gaussian_dynamics.analytic_soc_models_v220 import SingletTripletSOCConfigV220
from gaussian_dynamics.correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from gaussian_dynamics.curved_basis_adaptation_v293 import finite_difference_full_jacobian_v293, build_tdvp_snapshot_v293
from gaussian_dynamics.curved_geometry_v293 import straight_path_v293, wilson_path_v293
from gaussian_dynamics.molecular_frame_bridge_v300 import MolecularFrameBridgeV300
from gaussian_dynamics.molecular_hamiltonian_v301 import build_molecular_grid_hamiltonian_v301, build_molecular_tdvp_snapshot_v301


def pair():
    return CorrelatedGaussianSpinorStateV270(
        [[-.7],[.8]], [[.3],[-.2]], [[[1.3]],[[1.1]]], [[[.07]],[[-.02]]],
        [[.8,.1j,.05,.03j],[.05j,.2,.07j,.01]],
    ).normalized()


@pytest.mark.parametrize('moving', [False, True])
def test_exact_transport_matches_independent_nac_ode_and_frame_endpoints(moving):
    _,bridge,_,gauge,_=analytic_soc_fixture_v307(moving=moving)
    for start,end in (([-.71],[.92]),([.2],[-1.1]),([0.],[0.])):
        path=straight_path_v293(start,end)
        exact=wilson_path_v293(bridge.connection_field(),path)
        ode=wilson_path_v293(bridge.numerical_connection_field(),path)
        expected=gauge.matrix(end).conj().T@gauge.matrix(start) if moving else np.eye(4)
        np.testing.assert_allclose(exact,ode,atol=3e-12,rtol=0)
        np.testing.assert_allclose(exact,expected,atol=3e-13,rtol=0)


@pytest.mark.parametrize('moving', [False, True])
def test_grid_operator_and_tangent_agree_with_numerical_molecular_bridge(moving):
    state,bridge,grid,_,_=analytic_soc_fixture_v307(moving=moving)
    numerical=MolecularFrameBridgeV300(bridge.provider,[0.])
    reference=build_molecular_grid_hamiltonian_v301(numerical,grid.axes)
    a=build_molecular_tdvp_snapshot_v301(state,bridge,grid)
    b=build_molecular_tdvp_snapshot_v301(state,numerical,reference)
    np.testing.assert_allclose(grid.matrix.toarray(),reference.matrix.toarray(),atol=3e-12,rtol=0)
    np.testing.assert_allclose(a.snapshot.jacobian,b.snapshot.jacobian,atol=2e-9,rtol=0)
    np.testing.assert_allclose(a.snapshot.velocity,b.snapshot.velocity,atol=1e-8,rtol=0)
    assert grid.fingerprint != reference.fingerprint


@pytest.mark.parametrize('moving', [False, True])
def test_analytic_tangent_has_second_order_central_difference_agreement(moving):
    _,bridge,grid,gauge,_=analytic_soc_fixture_v307(moving=moving)
    state=pair()
    if moving:
        state=replace(state,coefficients=np.array([gauge.matrix(q).conj().T@c for q,c in zip(state.q,state.coefficients)]))
    connection=bridge.connection_field()
    exact=bridge.parameter_jacobian(state,connection,grid.quadrature,path_factory=straight_path_v293,active_shape_mask=None)
    errors=[]
    for h in (.004,.002,.001):
        finite=finite_difference_full_jacobian_v293(state,connection,grid.quadrature,relative_step=h)
        errors.append(np.linalg.norm(exact-finite)/np.linalg.norm(exact))
    assert errors[-1]<2e-6
    assert all(3.8<a/b<4.2 for a,b in zip(errors[:-1],errors[1:]))


def test_full_shape_velocity_and_physical_time_derivative_are_gauge_covariant():
    _,fixed,grid,gauge,_=analytic_soc_fixture_v307()
    _,moving,moving_grid,_,_=analytic_soc_fixture_v307(moving=True)
    state=pair()
    transformed=replace(state,coefficients=np.array([gauge.matrix(q).conj().T@c for q,c in zip(state.q,state.coefficients)]))
    a=build_molecular_tdvp_snapshot_v301(state,fixed,grid).snapshot
    b=build_molecular_tdvp_snapshot_v301(transformed,moving,moving_grid).snapshot
    np.testing.assert_allclose(a.velocity[16:],b.velocity[16:],atol=2e-9,rtol=0)
    da=(a.jacobian@a.velocity).reshape(-1,4)
    db=(b.jacobian@b.velocity).reshape(-1,4)
    db=np.array([gauge.matrix(q)@v for q,v in zip(grid.quadrature.points,db)])
    np.testing.assert_allclose(da,db,atol=2e-10,rtol=0)


def test_analytic_provider_configuration_mutation_fails_closed():
    state,bridge,grid,_,_=analytic_soc_fixture_v307()
    bridge.provider.config=replace(bridge.provider.config,mass=3.)
    with pytest.raises(ValueError,match='changed'):
        build_molecular_tdvp_snapshot_v301(state,bridge,grid)


def test_gauge_replacement_and_connection_corruption_fail_closed(monkeypatch):
    state,bridge,grid,gauge,_=analytic_soc_fixture_v307(moving=True)
    bridge.provider.gauge=gauge
    with pytest.raises(ValueError,match='changed'):
        bridge.connection_field()
    _,bridge,grid,_,_=analytic_soc_fixture_v307()
    original=bridge.provider.evaluate_snapshot
    def corrupt(q):
        snapshot=original(q)
        return replace(snapshot,point=replace(snapshot.point,connection_q=1j*np.eye(4)[None,:,:]))
    monkeypatch.setattr(bridge.provider,'evaluate_snapshot',corrupt)
    with pytest.raises(ValueError,match='constant analytic connection'):
        bridge.operator_point([.3])


@pytest.mark.parametrize('bad', [np.zeros((1,1)),np.full((52,12),np.nan)])
def test_supplied_tangent_must_be_complete_finite_and_weighted_shape(bad):
    state,bridge,grid,_,_=analytic_soc_fixture_v307()
    with pytest.raises(ValueError,match='Jacobian'):
        build_tdvp_snapshot_v293(state,bridge.connection_field(),grid.quadrature,np.zeros((13,4)),jacobian=bad)


@pytest.mark.parametrize('config,gauge', [(object(),None),(SingletTripletSOCConfigV220(),object())])
def test_fast_bridge_cannot_claim_an_arbitrary_provider(config,gauge):
    with pytest.raises(TypeError):
        AnalyticSOCFrameBridgeV307(config,gauge)
