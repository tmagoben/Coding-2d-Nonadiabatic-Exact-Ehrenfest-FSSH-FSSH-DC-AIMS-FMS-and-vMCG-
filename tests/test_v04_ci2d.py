import numpy as np

from gaussian_dynamics.ci2d import (
    LVC2DParameters,
    diabatic_potential_2d,
    adiabatic_energies_2d,
    adiabatic_gradients_2d,
    vector_nac_2d,
    circle_path,
    berry_line_integral,
    parallel_transport_real_state,
)


def test_ci_degeneracy_and_analytic_energies():
    p=LVC2DParameters()
    E0=adiabatic_energies_2d(0.0,0.0,p)
    assert abs(E0[1]-E0[0]) < 1e-15

    R=(0.7,-0.4)
    numeric=np.linalg.eigvalsh(diabatic_potential_2d(*R,p))
    analytic=adiabatic_energies_2d(*R,p)
    assert np.allclose(numeric,analytic,atol=1e-14)


def test_adiabatic_gradients_match_finite_difference():
    p=LVC2DParameters()
    R=np.array([0.8,0.5])
    g=adiabatic_gradients_2d(R,p)
    h=1e-6

    for a in (0,1):
        for k in (0,1):
            Rp=R.copy(); Rm=R.copy()
            Rp[k]+=h; Rm[k]-=h
            fd=(adiabatic_energies_2d(*Rp,p)[a]-adiabatic_energies_2d(*Rm,p)[a])/(2*h)
            assert abs(fd-g[a,k]) < 2e-9


def test_vector_nac_is_antisymmetric():
    d=vector_nac_2d((0.7,0.4))
    assert np.allclose(d[0,1],-d[1,0],atol=1e-14)


def test_berry_line_integral_is_pi_for_one_winding():
    path=circle_path(radius=1.3,n=4001)
    value=berry_line_integral(path)
    assert abs(value-np.pi) < 2e-6


def test_continuous_real_transport_returns_minus_sign():
    path=circle_path(radius=1.1,n=1001)
    _, overlap=parallel_transport_real_state(path,state=0)
    assert overlap < -0.999
