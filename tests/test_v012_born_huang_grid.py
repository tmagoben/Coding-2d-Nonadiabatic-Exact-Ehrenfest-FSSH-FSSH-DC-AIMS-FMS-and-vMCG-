import numpy as np

from gaussian_dynamics.benchmark_campaign import CIPassageConfig
from gaussian_dynamics.dynamic_graph_aims import DynamicGraphTBF
from gaussian_dynamics.exact_benchmark import localized_adiabatic_packet_2d
from gaussian_dynamics.born_huang_grid_v12 import (
    build_born_huang_grid_2d,
    build_born_huang_matrices,
    reconstruct_born_huang_wavefunction,
    born_huang_reduced_density,
    born_huang_basis_time_matrix_grid,
)


def _seed(config):
    return DynamicGraphTBF(
        uid=0,
        state=config.state,
        q=config.q_array(),
        p=config.p_array(),
        A=config.A_matrix(),
        node=("seed",0),
    )


def test_single_born_huang_tbf_reconstructs_exact_initial_packet_on_same_grid():
    config=CIPassageConfig()
    grid=build_born_huang_grid_2d(
        grid_n=40,
        half_width=config.half_width,
        mass=config.mass,
    )
    basis=[_seed(config)]

    psi=reconstruct_born_huang_wavefunction(
        np.array([1.0+0j]),
        basis,
        grid,
    )
    exact=localized_adiabatic_packet_2d(
        grid.points,
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        state=config.state,
    )

    assert np.allclose(psi,exact,atol=2e-13)


def test_born_huang_initial_reduced_density_matches_direct_grid_density():
    config=CIPassageConfig()
    grid=build_born_huang_grid_2d(
        grid_n=40,
        half_width=config.half_width,
        mass=config.mass,
    )
    basis=[_seed(config)]

    rho=born_huang_reduced_density(
        [1.0+0j],basis,grid,normalize=True
    )

    psi=localized_adiabatic_packet_2d(
        grid.points,
        config.q_array(),
        config.p_array(),
        config.A_matrix(),
        state=config.state,
    )
    flat=psi.reshape(-1,2)
    direct=(flat.T@np.conj(flat))*grid.area
    direct/=np.trace(direct)

    assert np.allclose(rho,direct,atol=2e-13)


def test_projected_born_huang_S_H_are_hermitian():
    config=CIPassageConfig()
    grid=build_born_huang_grid_2d(
        grid_n=28,
        half_width=config.half_width,
        mass=config.mass,
    )

    basis=[
        _seed(config),
        DynamicGraphTBF(
            uid=1,
            state=0,
            q=np.array([-0.25,0.35]),
            p=np.array([8.0,0.5]),
            A=np.diag([0.9,1.4]),
            node=("x",1),
        ),
    ]

    S,H=build_born_huang_matrices(basis,grid)

    assert np.allclose(S,S.conj().T,atol=1e-12)
    assert np.allclose(H,H.conj().T,atol=1e-12)
    assert np.min(np.linalg.eigvalsh(S)) > 0.0


def test_grid_basis_time_matrix_matches_finite_difference_of_grid_overlap():
    config=CIPassageConfig()
    grid=build_born_huang_grid_2d(
        grid_n=30,
        half_width=config.half_width,
        mass=config.mass,
    )

    basis=[
        DynamicGraphTBF(
            uid=0,state=1,
            q=np.array([0.6,0.5]),
            p=np.array([0.4,0.2]),
            A=np.diag([1.2,0.9]),
            node=("a",0),
        ),
        DynamicGraphTBF(
            uid=1,state=1,
            q=np.array([0.2,0.8]),
            p=np.array([-0.1,0.35]),
            A=np.diag([0.8,1.3]),
            node=("b",1),
        ),
    ]

    qdots=np.array([[0.07,-0.03],[-0.02,0.05]])
    pdots=np.array([[0.01,0.02],[-0.03,0.01]])

    T=born_huang_basis_time_matrix_grid(
        basis,grid,qdots,pdots
    )

    h=1e-6

    def shifted(sign):
        out=[]
        for i,b in enumerate(basis):
            out.append(
                DynamicGraphTBF(
                    uid=b.uid,state=b.state,
                    q=b.q+sign*h*qdots[i],
                    p=b.p+sign*h*pdots[i],
                    A=b.A.copy(),
                    node=b.node,
                )
            )
        return out

    Sp,_=build_born_huang_matrices(shifted(+1),grid)
    Sm,_=build_born_huang_matrices(shifted(-1),grid)
    Sdot=(Sp-Sm)/(2*h)

    assert np.allclose(Sdot,T+T.conj().T,atol=5e-7)
