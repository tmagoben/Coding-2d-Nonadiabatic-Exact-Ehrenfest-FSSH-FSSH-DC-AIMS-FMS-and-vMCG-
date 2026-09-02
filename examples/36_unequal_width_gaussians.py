import numpy as np

from gaussian_dynamics.gaussian_general import (
    gaussian_overlap_general,
    gaussian_cross_centroid,
    real_overlap_saddle_point,
    kinetic_matrix_element_general,
)

qi=np.array([-0.6,0.3])
pi=np.array([0.7,-0.15])
Ai=np.array([[1.25,0.10],[0.10,0.80]])

qj=np.array([0.5,-0.45])
pj=np.array([-0.2,0.65])
Aj=np.array([[0.75,-0.08],[-0.08,1.10]])

M=np.array([[2.7,0.2],[0.2,1.9]])

print("Unequal-width Gaussian algebra")
print("------------------------------")
print("<g_i|g_j> =",gaussian_overlap_general(qi,pi,Ai,qj,pj,Aj))
print("real overlap saddle =",real_overlap_saddle_point(qi,Ai,qj,Aj))
print("complex cross centroid =",gaussian_cross_centroid(qi,pi,Ai,qj,pj,Aj))
print("<g_i|T|g_j> =",kinetic_matrix_element_general(
    qi,pi,Ai,qj,pj,Aj,M
))
