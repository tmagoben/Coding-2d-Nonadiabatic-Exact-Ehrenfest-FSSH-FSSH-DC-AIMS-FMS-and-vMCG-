import numpy as np

from gaussian_dynamics.ci2d import analytic_adiabatic_vectors
from gaussian_dynamics.gauge_graph import ElectronicGaugeGraph

n=80
angles=np.linspace(-np.pi,np.pi,n,endpoint=False)
states=[]
for angle in angles:
    R=np.array([np.cos(angle),np.sin(angle)])
    states.append(analytic_adiabatic_vectors(R)[:,0])

graph=ElectronicGaugeGraph(1)
for i in range(n):
    j=(i+1)%n
    O=np.array([[np.vdot(states[i],states[j])]])
    graph.add_overlap(i,j,O)

W=graph.wilson_loop(list(range(n)))

tree_gauges,_=graph.spanning_tree_gauges(0)
sync_gauges=graph.synchronize(0,max_iter=100,tolerance=1e-10,restarts=2)

print("CI lower-state gauge graph")
print("--------------------------")
print("number of nodes:",n)
print("Wilson loop:",W[0,0])
print("Berry phase / pi:",np.angle(W[0,0])/np.pi)
print("tree-gauge objective:",graph.gauge_objective(tree_gauges))
print("synchronized objective:",graph.gauge_objective(sync_gauges))
print("\nThe nontrivial W=-1 holonomy remains after gauge smoothing.")
