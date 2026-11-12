"""Exact constant-connection transport for the declared analytic SOC fixture.

This is deliberately restricted to the complete singlet/triplet model and the
linear-phase gauge. It is not an acceleration assumption for arbitrary molecular
providers. Independent NAC-ODE comparisons remain part of validation.
"""
from dataclasses import asdict, replace

import numpy as np
from scipy.linalg import expm

from .analytic_soc_models_v220 import AnalyticSingletTripletSOCProviderV220, SingletTripletSOCConfigV220
from .complex_gauge_v21 import GaugeTransformedOperatorProviderV21, PhaseMixingGaugeV21, random_unitary_v21
from .correlated_gaussian_tdvp_v270 import CorrelatedGaussianSpinorStateV270
from .curved_geometry_v293 import ConnectionFieldV293
from .curved_basis_adaptation_v293 import coefficient_basis_v293, active_parameter_indices_v293
from .molecular_frame_bridge_v300 import MolecularFrameBridgeV300, MolecularFrameBridgeSettingsV300, _hash_v300
from .molecular_hamiltonian_v301 import build_molecular_grid_hamiltonian_v301

SCHEMA_V307_TRANSPORT = "gnd-analytic-soc-constant-connection-v0.30.7-dev"
SCHEMA_V307_TANGENT = "gnd-analytic-constant-connection-1d-tangent-v0.30.7-dev"


class AnalyticSOCFrameBridgeV307(MolecularFrameBridgeV300):
    """Integrate A=0 or A=U0† diag(i B) U0 exactly with a matrix exponential."""
    def __init__(self, config, gauge=None, *, settings=MolecularFrameBridgeSettingsV300()):
        if type(config) is not SingletTripletSOCConfigV220:
            raise TypeError("exact analytic SOC bridge requires SingletTripletSOCConfigV220")
        config = replace(config).validate()
        if gauge is not None:
            if type(gauge) is not PhaseMixingGaugeV21 or (gauge.nq, gauge.nstate) != (1, 4):
                raise TypeError("exact analytic SOC bridge requires a four-state one-dimensional PhaseMixingGaugeV21")
            gauge = PhaseMixingGaugeV21(gauge.U0.copy(), gauge.phase_gradient.copy(), gauge.phase_offset.copy())
            for value in (gauge.U0, gauge.phase_gradient, gauge.phase_offset):
                value.setflags(write=False)
        provider = AnalyticSingletTripletSOCProviderV220(config)
        if gauge is not None:
            provider = GaugeTransformedOperatorProviderV21(provider, gauge)
        self.gauge = gauge
        self._generators = (np.zeros((1, 4, 4), complex) if gauge is None else
                            np.array([gauge.U0.conj().T @ np.diag(1j*gauge.phase_gradient[:, 0]) @ gauge.U0]))
        self._generators.setflags(write=False)
        super().__init__(provider, [0.0], settings=settings)
        self._analytic_identity = self._identity()
        self.bridge_fingerprint = _hash_v300({"schema": SCHEMA_V307_TRANSPORT,
                                              "parent_bridge": self.bridge_fingerprint,
                                              "tangent_backend": SCHEMA_V307_TANGENT,
                                              "analytic_identity": self._analytic_identity})

    def _identity(self):
        if self.gauge is None:
            base, gauge_payload = self.provider, None
        else:
            if type(self.provider) is not GaugeTransformedOperatorProviderV21 or self.provider.gauge is not self.gauge:
                raise ValueError("analytic SOC provider or gauge changed")
            base = self.provider.base_provider
            gauge_payload = {"U0": self.gauge.U0, "gradient": self.gauge.phase_gradient, "offset": self.gauge.phase_offset}
        if type(base) is not AnalyticSingletTripletSOCProviderV220:
            raise ValueError("analytic SOC provider changed")
        return _hash_v300({"config": asdict(base.config), "gauge": gauge_payload, "generators": self._generators})

    def _check_identity(self):
        if self._identity() != self._analytic_identity:
            raise ValueError("analytic SOC provider or gauge changed after bridge construction")

    def evaluate_snapshot(self, q):
        self._check_identity()
        snapshot = super().evaluate_snapshot(q)
        if not np.allclose(snapshot.point.connection_q, self._generators, rtol=0, atol=2e-12):
            raise ValueError("provider connection does not match the declared constant analytic connection")
        return snapshot

    def connection_field(self, *, transport_mode="nac_ode"):
        self._check_identity()
        if transport_mode != "nac_ode":
            return super().connection_field(transport_mode=transport_mode)
        def transport(start, end):
            self._check_identity()
            return expm(-np.einsum("m,mab->ab", np.asarray(end)-np.asarray(start), self._generators))
        return ConnectionFieldV293(1, 4, lambda q: self.operator_point(q).connection_q,
                                  SCHEMA_V307_TRANSPORT, transport).validate([self.reference_q])

    def numerical_connection_field(self):
        """Independent NAC-ODE route; no analytic segment callback."""
        return MolecularFrameBridgeV300.connection_field(self, transport_mode="nac_ode")

    def parameter_jacobian(self, state, connection, quadrature, *, path_factory, active_shape_mask):
        """Exact coefficient, q, p, log-width and chirp derivatives in 1D.

        W(R<-q)=exp[-A(R-q)] gives dW/dq=W A. The center derivative
        therefore includes a transported A c term in addition to dg/dq.
        """
        self._check_identity()
        state.validate(require_normalized=False)
        if state.ndim != 1 or state.nstate != 4:
            raise ValueError("analytic SOC tangent requires a one-dimensional four-state spinor")
        phi = coefficient_basis_v293(state, connection, quadrature, path_factory=path_factory)
        weight = np.repeat(np.sqrt(quadrature.weights), state.nstate)
        coefficients = phi.reshape(quadrature.npoint*state.nstate, -1)
        blocks = [coefficients, 1j*coefficients]
        shapes = [[] for _ in range(4)]
        for i in range(state.ngaussian):
            x = quadrature.points[:, 0] - state.q[i, 0]
            width, chirp, momentum = state.width_matrices[i, 0, 0], state.chirp_matrices[i, 0, 0], state.p[i, 0]
            field = np.einsum("nsa,a->ns", phi[:, :, i, :], state.coefficients[i])
            connection_term = np.einsum("nsa,a->ns", phi[:, :, i, :], self._generators[0] @ state.coefficients[i])
            derivatives = [(((width-1j*chirp)*x-1j*momentum)[:, None]*field + connection_term),
                           1j*x[:, None]*field, (.25-.5*width*x*x)[:, None]*field,
                           .5j*(x*x)[:, None]*field]
            for group, derivative in zip(shapes, derivatives):
                group.append(derivative.reshape(-1))
        blocks.extend(np.stack(group, axis=1) for group in shapes)
        full = weight[:, None] * np.concatenate(blocks, axis=1)
        return full[:, active_parameter_indices_v293(state, active_shape_mask)]


def analytic_soc_fixture_v307(*, moving=False, soc_enabled=True, q0=.1, p0=.3):
    config = SingletTripletSOCConfigV220(mass=2.0, soc_scale=.2 if soc_enabled else 0.0, soc_enabled=soc_enabled)
    gauge = PhaseMixingGaugeV21(random_unitary_v21(4, seed=303),
                              np.array([[.3], [-.2], [.1], [.25]]), np.zeros(4))
    bridge = AnalyticSOCFrameBridgeV307(config, gauge if moving else None)
    grid = build_molecular_grid_hamiltonian_v301(bridge, (np.linspace(-3, 3, 13),))
    state = CorrelatedGaussianSpinorStateV270([[q0]], [[p0]], [[[1.1]]], [[[.08]]], [[1., 0., 0., 0.]]).normalized()
    if moving:
        state = replace(state, coefficients=np.array([gauge.matrix(state.q[0]).conj().T @ state.coefficients[0]]))
    return state, bridge, grid, gauge, config
