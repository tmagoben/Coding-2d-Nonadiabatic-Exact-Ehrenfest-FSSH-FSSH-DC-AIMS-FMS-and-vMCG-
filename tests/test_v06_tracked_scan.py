import numpy as np

from gaussian_dynamics.molecular_backend import GeneralizedElectronicStructurePoint
from gaussian_dynamics.tracked_scan import (
    run_tracked_scan,
    TrackedScan1DProvider,
)


class SimpleSequentialProvider:
    def evaluate(self,q):
        x=float(np.asarray(q).reshape(-1)[0])
        E=np.array([x,1.0-x])
        G=np.array([[1.0],[-1.0]])
        D=np.zeros((2,2,1))
        D[0,1,0]=0.2+x
        D[1,0,0]=-(0.2+x)
        return GeneralizedElectronicStructurePoint(
            q=np.array([x]),
            energies=E,
            gradients_q=G,
            nac_q=D,
            mass_matrix_q_au=np.array([[10.0]]),
            metadata={"x":x},
        ).validate()


def test_tracked_scan_and_order_independent_interpolation():
    q=np.array([-1.0,0.0,1.0])
    scan=run_tracked_scan(q,SimpleSequentialProvider())

    assert scan.q.shape==(3,1)
    assert scan.energies.shape==(3,2)

    provider=TrackedScan1DProvider(scan)

    p1=provider.evaluate(np.array([0.25]))
    p2=provider.evaluate(np.array([-0.25]))
    p3=provider.evaluate(np.array([0.25]))

    assert np.allclose(p1.energies,p3.energies)
    assert np.allclose(p1.nac_q,p3.nac_q)
    assert np.isclose(p1.nac_q[0,1,0],0.45)
    assert np.isclose(p2.nac_q[0,1,0],-0.05)
