import numpy as np

from gaussian_dynamics.basis_management import prune_redundant_basis

# Nearly duplicate normalized basis functions.
eps=1e-10
S=np.array([
    [1.0,1.0-eps],
    [1.0-eps,1.0],
],dtype=complex)
C=np.array([0.5,0.5],dtype=complex)

out=prune_redundant_basis(
    C,S,
    condition_limit=1e8,
    eigenvalue_floor=1e-8,
    max_projection_loss=1e-6,
)

print('condition before:',out.condition_before)
print('condition after: ',out.condition_after)
print('removed indices: ',out.removed)
print('projection loss: ',out.projection_loss)
print('projected C:     ',out.coefficients)
