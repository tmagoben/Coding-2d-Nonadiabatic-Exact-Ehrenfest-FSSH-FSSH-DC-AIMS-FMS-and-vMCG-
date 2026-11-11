"""Audited covariant kinetic + electronic action for v0.30.1 development.

On a uniform tensor grid with zero Dirichlet ghost values, discretize

    H = -sum_a (partial_a + D_a)^2 / (2 m_a) + H_el(R).

Each nearest-neighbour link is a Wilson line of the provider's NAC connection,
checked against an independently certified overlap polar factor. Reverse blocks
are adjoints, so the discrete Hamiltonian is Hermitian and gauge covariant.
Only constant diagonal masses are admitted. This projected model does not infer
an omitted-state Born--Huang correction or certify a live molecular trajectory.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .curved_basis_adaptation_v292 import CurvedQuadratureV292, CurvedTDVPSnapshotV292
from .curved_basis_adaptation_v293 import (
    build_tdvp_snapshot_v293,
    evaluate_curved_state_v293,
)
from .curved_geometry_v293 import straight_path_v293, wilson_path_v293
from .molecular_frame_bridge_v300 import MolecularFrameBridgeV300, _hash_v300


V301_HAMILTONIAN_SCHEMA = "gnd-covariant-grid-hamiltonian-v0.30.1-dev"


def _uniform_axes(axes, ndim):
    result = []
    for axis in axes:
        raw = np.asarray(axis)
        if np.iscomplexobj(raw):
            raise ValueError("grid axes must be real")
        values = np.array(raw, dtype=float, copy=True)
        if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
            raise ValueError("grid axes must be finite vectors with at least two points")
        steps = np.diff(values)
        if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-11, atol=0.0):
            raise ValueError("grid axes must be strictly increasing and uniform")
        values.setflags(write=False)
        result.append(values)
    if len(result) != ndim:
        raise ValueError("grid and bridge coordinate dimensions disagree")
    return tuple(result)


@dataclass(frozen=True)
class MolecularGridHamiltonianV301:
    axes: tuple
    quadrature: CurvedQuadratureV292
    mass_matrix: np.ndarray
    electronic_blocks: np.ndarray
    kinetic: sparse.csr_matrix
    matrix: sparse.csr_matrix
    segment_audits: tuple
    point_fingerprints: tuple
    bridge_fingerprint: str
    fingerprint: str

    @property
    def nstate(self):
        return self.electronic_blocks.shape[1]

    def apply(self, psi):
        """Apply both nuclear kinetic energy and the total electronic operator."""
        psi = np.asarray(psi, dtype=complex)
        if psi.shape != (self.quadrature.npoint, self.nstate):
            raise ValueError("wavefunction must have shape (npoint, nstate)")
        if not np.all(np.isfinite(psi)):
            raise ValueError("wavefunction must be finite")
        return (self.matrix @ psi.reshape(-1)).reshape(psi.shape)

    def hamiltonian_action(self, state, psi, quadrature):
        """Adapter for the recovered curved-energy and basis-adaptation APIs."""
        state.validate(require_normalized=False)
        quadrature.validate()
        if state.ndim != len(self.axes) or state.nstate != self.nstate:
            raise ValueError("state and Hamiltonian dimensions disagree")
        if not (
            np.array_equal(quadrature.points, self.quadrature.points)
            and np.array_equal(quadrature.weights, self.quadrature.weights)
        ):
            raise ValueError("Hamiltonian action requires its own grid and weights")
        return self.apply(psi)

    def evidence(self):
        """Finite, serializable diagnostics; historical release status is separate."""
        defect = self.matrix - self.matrix.conj().T
        return {
            "schema": V301_HAMILTONIAN_SCHEMA,
            "fingerprint": self.fingerprint,
            "bridge_fingerprint": self.bridge_fingerprint,
            "grid_shape": [len(axis) for axis in self.axes],
            "nstate": self.nstate,
            "matrix_dimension": self.matrix.shape[0],
            "matrix_nnz": self.matrix.nnz,
            "boundary": "zero Dirichlet ghost values, one grid step beyond each axis",
            "mass_matrix_au": self.mass_matrix.tolist(),
            "audited_segments": len(self.segment_audits),
            "maximum_connection_overlap_residual": max(
                receipt.connection_overlap_residual for receipt in self.segment_audits
            ),
            "hermiticity_residual": float(sparse.linalg.norm(defect)),
            "point_fingerprints": list(self.point_fingerprints),
            "production_molecular_accuracy_certified": False,
        }


def build_molecular_grid_hamiltonian_v301(bridge, axes):
    """Sample a stationary bridge and assemble an audited sparse Hamiltonian.

    Axis entries are interior nodes; ghost values outside the box are zero.
    No periodic wrap link is introduced. The provider must describe rectilinear
    coordinates with a constant diagonal mass. Mass constancy is checked at the
    reference point and every grid point; it cannot be proved between samples.
    """
    if not isinstance(bridge, MolecularFrameBridgeV300):
        raise TypeError("bridge must be MolecularFrameBridgeV300")
    axes = _uniform_axes(axes, bridge.ndim)
    shape = tuple(len(axis) for axis in axes)
    spacing = np.array([axis[1] - axis[0] for axis in axes])
    mesh = np.meshgrid(*axes, indexing="ij")
    points = np.stack([part.reshape(-1) for part in mesh], axis=1)
    weights = np.full(len(points), np.prod(spacing))
    quadrature = CurvedQuadratureV292(points, weights).validate()

    reference_mass = bridge.operator_point(bridge.reference_q).mass_matrix_q_au.copy()
    diagonal = np.diag(np.diag(reference_mass))
    tolerance = 1e-12 * max(float(np.linalg.norm(reference_mass)), 1.0)
    if np.linalg.norm(reference_mass - diagonal) > tolerance:
        raise ValueError("v0.30.1 requires a diagonal mass matrix")
    electronic = []
    point_fingerprints = []
    for point in points:
        snapshot = bridge.operator_point(point)
        if np.linalg.norm(snapshot.mass_matrix_q_au - reference_mass) > tolerance:
            raise ValueError("v0.30.1 requires a constant mass matrix on the grid")
        electronic.append(snapshot.H.copy())
        point_fingerprints.append(_hash_v300({
            "provider_fingerprint": bridge.provider_fingerprint,
            "point": bridge._point_payload(snapshot),
        }))
    electronic = np.array(electronic)
    nstate = bridge.nstate
    npoint = len(points)
    prefactors = 1.0 / (2.0 * np.diag(reference_mass) * spacing**2)
    if not np.all(np.isfinite(prefactors)) or np.any(prefactors <= 0):
        raise ValueError("grid kinetic coefficients must be finite and positive")
    dimension = npoint * nstate
    rows = list(range(dimension))
    columns = list(range(dimension))
    data = [2.0 * float(np.sum(prefactors))] * dimension
    audits = []
    connection = bridge.connection_field(transport_mode="nac_ode")

    def add_block(left, right, block):
        r, c = np.indices((nstate, nstate))
        rows.extend((left * nstate + r).reshape(-1).tolist())
        columns.extend((right * nstate + c).reshape(-1).tolist())
        data.extend(block.reshape(-1).tolist())

    for index in np.ndindex(shape):
        left = int(np.ravel_multi_index(index, shape))
        for axis, prefactor in enumerate(prefactors):
            if index[axis] + 1 == shape[axis]:
                continue
            neighbour = list(index)
            neighbour[axis] += 1
            right = int(np.ravel_multi_index(tuple(neighbour), shape))
            audits.append(bridge.audit_segment(points[left], points[right]))
            # W maps the left frame into the right frame. In row left, col right
            # the coefficient vector must instead be transported right -> left.
            transport = wilson_path_v293(
                connection, straight_path_v293(points[left], points[right])
            )
            add_block(left, right, -prefactor * transport.conj().T)
            add_block(right, left, -prefactor * transport)

    kinetic = sparse.coo_matrix(
        (np.array(data, dtype=complex), (rows, columns)),
        shape=(dimension, dimension),
    ).tocsr()
    kinetic.eliminate_zeros()
    matrix = (kinetic + sparse.block_diag(electronic, format="csr")).tocsr()
    matrix.eliminate_zeros()
    fingerprint = _hash_v300({
        "schema": V301_HAMILTONIAN_SCHEMA,
        "bridge_fingerprint": bridge.bridge_fingerprint,
        "axes": axes,
        "boundary": "zero-dirichlet-ghosts",
        "mass": reference_mass,
        "point_fingerprints": point_fingerprints,
        "audits": [receipt.as_dict() for receipt in audits],
        "matrix_data": matrix.data,
        "matrix_indices": matrix.indices,
        "matrix_indptr": matrix.indptr,
    })
    for array in (quadrature.points, quadrature.weights, reference_mass, electronic):
        array.setflags(write=False)
    for array in (kinetic.data, kinetic.indices, kinetic.indptr,
                  matrix.data, matrix.indices, matrix.indptr):
        array.setflags(write=False)
    return MolecularGridHamiltonianV301(
        axes, quadrature, reference_mass, electronic, kinetic, matrix,
        tuple(audits), tuple(point_fingerprints), bridge.bridge_fingerprint,
        fingerprint,
    )


@dataclass(frozen=True)
class MolecularTDVPEvaluationV301:
    psi: np.ndarray
    hpsi: np.ndarray
    snapshot: CurvedTDVPSnapshotV292
    hamiltonian_fingerprint: str


def build_molecular_tdvp_snapshot_v301(
    state, bridge, hamiltonian, *, path_factory=straight_path_v293,
    relative_step=2e-6,
):
    """Evaluate the curved Gaussian state and feed its complete H psi to TDVP.

    The stationary provider must remain unchanged after grid assembly. Stored
    point fingerprints are rechecked before evaluating the tangent basis.
    The finite-difference tangent solve and its null-space gates are inherited
    from v0.29.3. This evaluates one state, not a time-propagation algorithm.
    """
    if bridge.bridge_fingerprint != hamiltonian.bridge_fingerprint:
        raise ValueError("Hamiltonian and molecular bridge fingerprints disagree")
    if not np.isfinite(relative_step) or relative_step <= 0:
        raise ValueError("relative_step must be finite and positive")
    for point, expected in zip(hamiltonian.quadrature.points, hamiltonian.point_fingerprints):
        observed = _hash_v300({
            "provider_fingerprint": bridge.provider_fingerprint,
            "point": bridge._point_payload(bridge.operator_point(point)),
        })
        if observed != expected:
            raise ValueError("molecular provider changed after Hamiltonian assembly")
    connection = bridge.connection_field(transport_mode="nac_ode")
    quadrature = hamiltonian.quadrature
    psi = evaluate_curved_state_v293(
        state, connection, quadrature, path_factory=path_factory,
    )
    hpsi = hamiltonian.apply(psi)
    snapshot = build_tdvp_snapshot_v293(
        state, connection, quadrature, hpsi,
        path_factory=path_factory, relative_step=relative_step,
    )
    return MolecularTDVPEvaluationV301(psi, hpsi, snapshot, hamiltonian.fingerprint)
