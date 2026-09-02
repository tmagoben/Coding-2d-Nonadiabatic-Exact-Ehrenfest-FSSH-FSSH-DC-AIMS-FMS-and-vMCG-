import numpy as np

from .ci2d import (
    LVC2DParameters,
    adiabatic_energies_2d,
    adiabatic_gradients_2d,
    vector_nac_2d,
)
from .molecular_backend import GeneralizedElectronicStructurePoint


class LVC2DGeneralizedProvider:
    """Backend-independent two-coordinate provider for v0.5 regression tests."""

    def __init__(self, nuclear_mass_au=20.0, params=LVC2DParameters()):
        self.nuclear_mass_au = float(nuclear_mass_au)
        self.params = params

    def evaluate(self, q):
        q = np.asarray(q, dtype=float)
        if q.shape != (2,):
            raise ValueError("LVC2DGeneralizedProvider expects q=(x,y).")

        E = adiabatic_energies_2d(q[0], q[1], self.params)
        grad = adiabatic_gradients_2d(q, self.params)
        nac = vector_nac_2d(q, self.params)

        return GeneralizedElectronicStructurePoint(
            q=q.copy(),
            energies=np.asarray(E, dtype=float),
            gradients_q=np.asarray(grad, dtype=float),
            nac_q=np.asarray(nac, dtype=float),
            mass_matrix_q_au=self.nuclear_mass_au*np.eye(2),
            metadata={
                "provider": "analytic_2d_lvc",
                "mass_au": self.nuclear_mass_au,
            },
        ).validate()
