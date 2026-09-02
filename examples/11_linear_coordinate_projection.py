import numpy as np

from gaussian_dynamics.molecular_backend import (
    MolecularGeometry,
    CartesianElectronicStructurePoint,
    LinearGeometryMap,
    GeneralizedCoordinateProvider,
)


class SimpleBackend:
    def evaluate(self, geometry):
        # Synthetic backend used only to demonstrate the projection algebra.
        E=np.array([0.0,0.1])

        grad=np.zeros((2,2,3))
        grad[0,1,2]=0.03
        grad[1,1,2]=-0.02

        nac=np.zeros((2,2,2,3))
        nac[0,1,1,2]=0.4
        nac[1,0,1,2]=-0.4

        return CartesianElectronicStructurePoint(
            geometry=geometry,
            energies=E,
            gradients_cart=grad,
            nac_cart=nac,
            masses_amu=np.array([7.0,1.0]),
            metadata={"example":"synthetic projection"},
        ).validate()


R0=np.array([
    [0.0,0.0,0.0],
    [0.0,0.0,3.0],
])

# q increases only the H z-coordinate in this transparent demonstration.
modes=np.zeros((1,2,3))
modes[0,1,2]=1.0

geomap=LinearGeometryMap(("Li","H"),R0,modes)
provider=GeneralizedCoordinateProvider(SimpleBackend(),geomap)

point=provider.evaluate(np.array([0.2]))

print("Projected generalized-coordinate point")
print("energies:",point.energies)
print("gradients_q:",point.gradients_q[:,0])
print("d_01/dq:",point.nac_q[0,1,0])
print("mass matrix (electron masses):",point.mass_matrix_q_au)
