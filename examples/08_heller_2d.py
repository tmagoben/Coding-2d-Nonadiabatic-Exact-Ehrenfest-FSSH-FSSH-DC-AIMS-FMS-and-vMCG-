import numpy as np

from gaussian_dynamics.heller_nd import run_thawed_gaussian_nd

mass=2.0
omega=0.7
q0=np.array([-0.8,0.5])
p0=np.array([0.4,-0.3])
sigma=np.sqrt(1.0/(2.0*mass*omega))

V=lambda q: 0.5*mass*omega**2*np.dot(q,q)
grad=lambda q: mass*omega**2*np.asarray(q)
hess=lambda q: mass*omega**2*np.eye(2)

out=run_thawed_gaussian_nd(
    q0,p0,sigma,mass,V,grad,hess,
    dt=0.001,steps=1000,store_every=100,
)

print("2D Heller TGA")
print("final q:",out["q"][-1])
print("final p:",out["p"][-1])
print("final width matrix A:\n",out["A"][-1])
print("maximum symmetry error:",
      np.max(np.abs(out["A"]-np.swapaxes(out["A"],-1,-2))))
