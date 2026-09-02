import numpy as np

from gaussian_dynamics.state_tracking import (
    maximum_overlap_assignment,
    transform_state_properties,
)

# Previous tracked states versus current raw roots.
# State 0 has become raw root 1 and also acquired a minus sign.
O=np.array([
    [0.05,-0.96],
    [0.94, 0.03],
])

tracking=maximum_overlap_assignment(
    O,
    minimum_overlap=0.5,
    minimum_score_margin=0.1,
)

raw_E=np.array([0.10,0.90])
raw_grad=np.array([[1.0],[2.0]])

raw_d=np.array([
    [0.0, 0.3],
    [-0.3,0.0],
])

E,G,d=transform_state_properties(
    raw_E,
    raw_grad,
    raw_d,
    tracking,
)

print("Raw overlap matrix:")
print(O)
print("\nTracked -> raw permutation:",tracking.permutation)
print("Phase/sign corrections:",tracking.phase_factors)
print("Assigned positive overlaps:",tracking.assigned_overlaps)

print("\nRaw energy order:",raw_E)
print("Tracked identity order:",E)

print("\nRaw NAC:")
print(raw_d)
print("Tracked/gauge-corrected NAC:")
print(d)
