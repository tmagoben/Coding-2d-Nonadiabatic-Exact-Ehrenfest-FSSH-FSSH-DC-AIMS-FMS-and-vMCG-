import numpy as np

from gaussian_dynamics.ci2d import diabatic_potential_2d
from gaussian_dynamics.gaussian_nd import gaussian_nd
from gaussian_dynamics.exact2d import run_exact_2d

n=48
L=5.0
dx=2*L/n
x=-L+(np.arange(n)+0.5)*dx
X,Y=np.meshgrid(x,x,indexing="ij")
points=np.stack([X,Y],axis=-1)

g=gaussian_nd(
    points,
    q=np.array([-2.0,0.7]),
    p=np.array([1.0,-0.1]),
    A=np.eye(2),
)

psi0=np.zeros((n,n,2),complex)
psi0[...,1]=g

V=diabatic_potential_2d(X,Y)

out=run_exact_2d(
    psi0,dx,dx,V,
    mass=20.0,dt=0.004,steps=200,store_every=20
)

print("Exact 2D two-state CI benchmark")
print("initial norm:",out["norm"][0])
print("final norm:",out["norm"][-1])
print("max norm drift:",np.max(np.abs(out["norm"]-1.0)))
print("final diabatic populations:",out["populations"][-1])
