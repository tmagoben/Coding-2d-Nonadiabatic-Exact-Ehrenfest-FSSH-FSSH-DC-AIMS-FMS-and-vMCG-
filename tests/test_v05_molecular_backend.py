import numpy as np

from gaussian_dynamics.molecular_backend import (
    MolecularGeometry,
    CartesianElectronicStructurePoint,
    LinearGeometryMap,
    GeneralizedCoordinateProvider,
    AMU_TO_ELECTRON_MASS,
)


class FakeCartesianBackend:
    def evaluate(self, geometry):
        # Two states, two atoms.
        energies=np.array([0.1,0.2])
        gradients=np.array([
            [[1.0,2.0,3.0],[-1.0,0.5,4.0]],
            [[-0.2,0.7,1.1],[0.4,-0.5,0.9]],
        ])

        nac=np.zeros((2,2,2,3))
        d=np.array([[0.1,-0.2,0.3],[-0.4,0.5,0.2]])
        nac[0,1]=d
        nac[1,0]=-d

        return CartesianElectronicStructurePoint(
            geometry=geometry,
            energies=energies,
            gradients_cart=gradients,
            nac_cart=nac,
            masses_amu=np.array([2.0,3.0]),
            metadata={"fake":True},
        ).validate()


def test_linear_coordinate_projection_and_mass_matrix():
    symbols=("A","B")
    R0=np.zeros((2,3))

    # Two independent generalized directions.
    modes=np.array([
        [[1.0,0.0,0.0],[0.0,0.5,0.0]],
        [[0.0,0.0,0.25],[0.0,0.0,1.0]],
    ])

    geomap=LinearGeometryMap(symbols,R0,modes)
    provider=GeneralizedCoordinateProvider(FakeCartesianBackend(),geomap)

    q=np.array([0.3,-0.2])
    point=provider.evaluate(q)

    raw=FakeCartesianBackend().evaluate(geomap.geometry(q))
    J=geomap.J

    expected_grad=raw.gradients_cart.reshape(2,-1) @ J
    expected_nac=np.einsum("ijr,ra->ija",raw.nac_cart.reshape(2,2,-1),J)

    masses_cart=np.repeat(np.array([2.0,3.0])*AMU_TO_ELECTRON_MASS,3)
    expected_M=J.T @ np.diag(masses_cart) @ J

    assert np.allclose(point.gradients_q,expected_grad)
    assert np.allclose(point.nac_q,expected_nac)
    assert np.allclose(point.mass_matrix_q_au,expected_M)
    assert np.min(np.linalg.eigvalsh(point.mass_matrix_q_au)) > 0
