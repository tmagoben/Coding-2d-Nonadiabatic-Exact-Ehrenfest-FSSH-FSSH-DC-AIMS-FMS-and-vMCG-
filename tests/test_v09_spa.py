import numpy as np

from gaussian_dynamics.gaussian_nd import gaussian_nd
from gaussian_dynamics.spa_matrix_elements import scalar_spa_matrix_element
from gaussian_dynamics.dynamic_gauge_graph import AnalyticCI2DFrameProvider, IncrementalElectronicGraph
from gaussian_dynamics.graph_gaussian import GraphGaussianTBF
from gaussian_dynamics.spa_matrix_elements import build_graph_gaussian_matrices_spa


def test_spa1_is_exact_for_linear_scalar_function():
    n=150; L=7.0; dx=2*L/n
    x=-L+(np.arange(n)+0.5)*dx
    X,Y=np.meshgrid(x,x,indexing='ij')
    P=np.stack([X,Y],axis=-1)

    A=np.array([[1.1,0.1],[0.1,0.9]])
    qi=np.array([-0.7,0.3]); pi=np.array([0.6,-0.2])
    qj=np.array([0.4,-0.5]); pj=np.array([-0.1,0.7])
    gi=gaussian_nd(P,qi,pi,A); gj=gaussian_nd(P,qj,pj,A)

    f0=0.37; grad=np.array([0.21,-0.14])
    qc=0.5*(qi+qj)
    f=f0+grad[0]*(X-qc[0])+grad[1]*(Y-qc[1])
    numeric=np.vdot(gi,f*gj)*dx*dx

    spa0=scalar_spa_matrix_element(qi,pi,qj,pj,A,f0,order=0)
    spa1=scalar_spa_matrix_element(qi,pi,qj,pj,A,f0,gradient_at_saddle=grad,order=1)

    assert abs(spa1-numeric) < 2e-10
    assert abs(spa0-numeric) > 1e-5


def test_graph_spa1_matrix_is_hermitian():
    provider=AnalyticCI2DFrameProvider()
    manager=IncrementalElectronicGraph(2)
    q1=np.array([0.7,0.5]); q2=np.array([1.0,0.6]); qc=0.5*(q1+q2)
    manager.add_from_provider('a',q1,provider)
    manager.add_from_provider('b',q2,provider,connect_to=['a'])
    manager.add_from_provider('c',qc,provider,connect_to=['a','b'])

    A=np.eye(2)
    e0=np.array([1,0],complex); e1=np.array([0,1],complex)
    basis=[
        GraphGaussianTBF('a',q1,np.array([0.3,0.1]),A,e0),
        GraphGaussianTBF('b',q2,np.array([-0.2,0.2]),A,e1),
    ]
    refs={(0,0):'a',(1,1):'b',(0,1):'c',(1,0):'c'}
    M=20*np.eye(2)
    S,H=build_graph_gaussian_matrices_spa(
        basis,manager.registry,M,lambda i,j:refs[(i,j)],order=1
    )
    assert np.allclose(S,S.conj().T,atol=1e-12)
    assert np.allclose(H,H.conj().T,atol=1e-11)
