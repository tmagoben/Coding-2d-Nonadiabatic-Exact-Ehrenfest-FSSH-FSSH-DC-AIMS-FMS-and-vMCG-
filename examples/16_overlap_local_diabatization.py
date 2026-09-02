import numpy as np

from gaussian_dynamics.overlap_transport import (
    current_to_previous_procrustes,
    directional_nac_from_overlap,
)

theta=0.25

# O_ij = <previous_i | current_j>
O=np.array([
    [np.cos(theta),-np.sin(theta)],
    [np.sin(theta), np.cos(theta)],
],dtype=complex)

Q,aligned,s=current_to_previous_procrustes(O)

print("Electronic overlap matrix:")
print(O)

print("\nCurrent-basis Procrustes rotation Q:")
print(Q)

print("\nOverlap after local diabatic alignment:")
print(aligned)

print("\nPrincipal-overlap singular values:")
print(s)

ds=1e-4
k=0.7
th=k*ds
Osmall=np.array([
    [np.cos(th),-np.sin(th)],
    [np.sin(th), np.cos(th)],
])

d=directional_nac_from_overlap(Osmall,ds)

print("\nDirectional NAC recovered from a small overlap step:")
print(d)
print("Expected off-diagonal magnitude:",k)
